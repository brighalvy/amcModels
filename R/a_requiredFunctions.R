# calculates log(1 - exp(x))
log1mexp <- function(x) {
  stopifnot(x < 0.0)
  if (x > -0.693147) {
    # .693147 ~= log(2.0))
    out <- log(-expm1(x))
  } else {
    out <- log1p(-exp(x))
  }
  out
}

## Map log(alpha) to log(z) using stick breaking methodology (Neal 2003)
## input a is on the log scale
alpha_to_z <- function(a) {
  L <- length(a) - 1
  z <- c()
  #z_1m <- c()
  for (l in 1:L) {
    if (l == 1) {
      z[l] <- a[l]
    } else{
      z[l] <- a[l] - log(1 - sum(exp(a[1:(l - 1)])))
    }
    if (z[l] > 0) {
      z[l] <- -1e-300
    }
    #z_1m[l] <- log1mexp(z[l])
  }
  return(z)
}

## Update draws of z using a pseudo-slice sampler:
## z is a vector of log(z)
update_z <- function(z, z_1m, J, n_i, K, ga, prior.alpha) {
  C <- length(z)
  if (sum(is.na(n_i[1, ])) == 0) {
    x <- (apply(n_i, 2, sum) + 1) / sqrt(sum((apply(n_i, 2, sum) + 1) ^ 2)) *
      4
  } else{
    na.cols <- which(is.na(n_i[1, ]))
    x <- (apply(n_i[, -na.cols], 2, sum) + 1) / sqrt(sum((apply(n_i[, -na.cols], 2, sum) +
                                                            1) ^ 2)) * 4
  }
  z_new <- c()
  z_1m_new <- c()
  a <- x[1:C]
  b <- rev(cumsum(rev(x)))[-1]
  psi <- pbeta(exp(z), a, b)
  # New proposal distribution:
  y <- l.alpha.f.cond(z, z_1m, J, n_i, K, ga, prior.alpha) - sum((dbeta(exp(z), a, b, log = TRUE))) + log(runif(1))
  # draw value of psi:
  L <- rep(0, C)
  R <- rep(1, C)
  psi_new <- c()
  psi_new <- runif(C)
  z_new <- log(qbeta(psi_new, a, b))
  z_1m_new <- sapply(z_new, log1mexp)
  while (y > (l.alpha.f.cond(z_new, z_1m_new, J, n_i, K, ga, prior.alpha) - sum((dbeta(
    exp(z_new), a, b, log = TRUE
  ))))) {
    L[psi_new < psi] <- psi_new[psi_new < psi]
    R[psi_new > psi] <- psi_new[psi_new > psi]
    psi_new <- runif(C, L, R)
    z_new <- log(qbeta(psi_new, a, b))
    z_1m_new <- sapply(z_new, log1mexp)
  }
  return(z_new)
}

## map z values to a:
## All done on the log scale
alpha_map <- function(z, z_1m) {
  L <- length(z) + 1
  a <- c()
  for (l in 1:L) {
    if (l == 1) {
      a[l] <- z[l]
    } else if (l != L) {
      a[l] <- z[l] + sum(z_1m[1:(l - 1)])
    } else {
      a[l] <- sum(z_1m)#log(1 - sum(exp(a[1:(l-1)])))
    }
  }
  return(a)
}

# Gamma slice sampler:
gamma_update <- function(z, z_1m, J, n_i, k, ga, g.a, g.b) {
  ## CDF Transformation:
  y <-  l.gamma.f.cond(z, z_1m, J, n_i, k, ga, g.a, g.b) + log(runif(1)) - pgamma(exp(ga), g.a, g.b)
  L <- 0
  R <- 1
  u <- runif(1, L, R)
  g <- log(qgamma(u, g.a, g.b))
  while (y >= l.gamma.f.cond(z, z_1m, J, n_i, k, g, g.a, g.b) - pgamma(exp(ga), g.a, g.b)) {
    if (g < ga) {
      L = u
    } else{
      R = u
    }
    u <- runif(1, L, R)
    g <- log(qgamma(u, g.a, g.b))
  }
  return(g)
}

# Log Full-Conditional for alpha:
## n is the data for the ith row with dimensions (K, j)
## all inputs (except n and J/K) are log values
l.alpha.f.cond <- function(z, z_1m, J, n_i, K, ga, prior.alpha) {
  k_obj <- c()
  alpha <- alpha_map(z, z_1m)
  for (k in 1:K) {
    k_obj[k] <- sum(lgamma(n_i[k, ][!is.na(n_i[k, ])] + exp(ga + alpha)) - lgamma(exp(ga +
                                                                                        alpha)))
  }
  # uses z ~ beta priors
  ## Set up priors:
  b <- c()
  if (prior.alpha == 0) {
    for (j in 1:(J - 1)) {
      b[j] <- (1 / J) * (J - j)
    }
    a <- 1 / J
  } else {
    for (j in 1:(J - 1)) {
      b[j] <- prior.alpha * (J - j)
    }
    a <- prior.alpha
  }
  res <- sum((a - 1) * (z)) + sum((b - 1) * (z_1m)) + sum(k_obj) - log(100)
  return(res)
}

#Log Full-Conditional for gamma:
l.gamma.f.cond <- function(z, z_1m, J, n_i, K, ga, g.a, g.b) {
  k_obj <- c()
  alpha <- alpha_map(z, z_1m)
  for (k in 1:K) {
    k_obj[k] <- lgamma(exp(ga)) + sum(lgamma(n_i[k, ][!is.na(n_i[k, ])] + exp(ga +
                                                                                alpha))) - lgamma(sum(n_i[k, ], na.rm = T) + exp(ga)) - sum(lgamma(exp(ga +
                                                                                                                                                         alpha)))
  }
  ## gamma priors are g.a, g.b
  res <- (g.a - 1) * ga - exp(ga) * g.b + sum(k_obj)
  return(res)
}

## Log EPA prior:
log_epa_prior <- function(p, alpha, delta, dist, sigma) {
  if (alpha < -delta | delta < 0 | delta >= 1) {
    # Set constraints
    -Inf
  }
  else {
    p_tm1 <- p[sigma[1]]
    q_tm1 <- 1
    count <- 2
    log_p_t <- log(1)
    for (i in sigma[2:length(p)]) {
      if (p[i] %in% p_tm1) {
        log_p_t[count] <- log(length(p_tm1) - delta * q_tm1) - log(alpha + length(p_tm1)) +
          log(sum(dist[sigma[which(p_tm1 == p[i])], i])) - log(sum(dist[sigma[1:length(p_tm1)], i]))
      } else{
        log_p_t[count] <- log(alpha + delta * q_tm1) - log(alpha + length(p_tm1))
      }
      count <- count + 1
      p_tm1 <- c(p_tm1, p[i])
      q_tm1 <- length(unique(p_tm1))
    }
    sum(log_p_t)
  }
}

log_delta_prior <- function(delta) {
  log(.4 * dbeta(delta, .05, 20) + .6 * dbeta(delta, 2, 2.5))
}
# Combine counts function:
combine_counts <- function(grouping, n_i) {
  ids <- (unique(grouping))
  res <- lapply(ids, function(id)
    sort(c(1:length(grouping))[grouping == id]))
  names(res) <- 1:length(res)
  counts.list <- lapply(res, \(x) {
    if (length(x) == 1) {
      n_i[x, ]
    }
    else{
      colSums(n_i[x, ])
    }
  })
  matrix(unlist(counts.list),
         ncol = ncol(n_i),
         byrow = T)
}
# Write full joint of gamma, alpha and counts
log_full_joint <- function(n, alpha, gamma, prior.alpha, g.a, g.b) {
  a <- matrix(alpha, nrow = 1)
  pa <- matrix(prior.alpha, nrow = 1)
  ag <- alpha * gamma
  if (is.null(dim(n))) {
    n <- matrix(n, nrow = 1)
  }
  n_ag <- t(apply(n, 1, \(x) x + ag))
  LaplacesDemon::ddirichlet((a), (pa), log = T) + (g.a - 1) * log(gamma) - g.b *
    gamma++nrow(n) * lgamma(gamma) + sum(-lgamma(rowSums(n_ag)) +
                                           rowSums(t(apply(
                                             lgamma(n_ag), 1, \(x) {
                                               x - lgamma(ag)
                                             }
                                           ))))
}
# Function to get posterior prob with given partition:
log_like_prob_group <- function(n_i,
                                p,
                                alpha,
                                gamma,
                                prior.alpha,
                                g.a,
                                g.b,
                                K) {
  ids <- (unique(p))
  res <- lapply(ids, function(id)
    sort(c(1:K)[p == id]))
  names(res) <- 1:length(res)
  counts.list <- lapply(res, \(x) {
    if (length(x) == 1) {
      n_i[x, ]
    }
    else{
      colSums(n_i[x, ])
    }
  })
  n_curr <- matrix(unlist(counts.list),
                   ncol = ncol(n_i),
                   byrow = T)
  log_full_joint(n_curr, alpha, gamma, prior.alpha, g.a, g.b)
}

# Function to update groupings using the all at once method:
update_groupings_aao <- function(p,
                                 beta,
                                 delta,
                                 dist,
                                 sigma,
                                 alpha,
                                 gamma,
                                 n_i,
                                 prior.alpha,
                                 g.a,
                                 g.b,
                                 K) {
  log.prior.probs <- apply(p, 1, \(x) log_epa_prior(x, beta, delta, dist, sigma))
  log.prior.probs <- log.prior.probs - max(log.prior.probs)
  ## Update groupings:
  log_probs <- apply(p,
                     1,
                     \(x) log_like_prob_group(n_i, x, alpha, gamma, prior.alpha, g.a, g.b, K))
  log_probs <- log_probs - max(log_probs)
  log_post_probs <- log_probs + log.prior.probs
  log_post_probs <- log_post_probs - max(log_post_probs)
  probs <- exp(log_post_probs) / sum(exp(log_post_probs))
  p[sample(1:52, 1, prob = probs), ]
}

# Function for sequential group update:
update_groupings_seq <- function(n_i,
                                 groupings,
                                 alpha,
                                 gamma,
                                 beta,
                                 delta,
                                 sigma,
                                 dist,
                                 prior.alpha,
                                 g.a,
                                 g.b,
                                 K) {
  count = 0
  for (i in sample(1:K)) {
    alloc <- groupings[-i]
    # Combine n's:
    ids <- (unique(alloc))
    res <- lapply(ids, function(id)
      sort(c(1:K)[-i][alloc == id]))
    names(res) <- 1:length(res)
    counts.list <- lapply(res, \(x) {
      if (length(x) == 1) {
        n_i[x, ]
      }
      else{
        colSums(n_i[x, ])
      }
    })
    n_curr <- matrix(unlist(counts.list),
                     nrow = length(res),
                     byrow = T)
    # Get log_post probabilities.
    log_probs <- log_lik <- log_prior <- c()
    for (j in 1:(length(ids) + 1)) {
      n_use <- n_curr
      p_use <- numeric(K)
      p_use[i] <- j
      for (g in 1:length(ids)) {
        p_use[res[[g]]] <- g
      }
      if (j != length(ids) + 1) {
        n_use[j, ] <- n_use[j, ] + n_i[i, ]
      } else{
        n_use <- rbind(n_use, n_i[i, ])
      }
      log_lik[j] <- log_full_joint(n_use, alpha, gamma, prior.alpha, g.a, g.b)
      log_prior[j] <- log_epa_prior(p_use, beta, delta, dist, sigma)
    }
    log_lik <- log_lik - max(log_lik)
    log_prior <- log_prior - max(log_prior)
    log_probs <- log_lik + log_prior
    log_probs <- log_probs - max(log_probs)
    probs <- exp(log_probs) / sum(exp(log_probs))
    n_g <- sample(1:(length(ids) + 1), 1, prob = probs)
    alloc.new <- numeric(K)
    for (j in 1:length(ids)) {
      alloc.new[res[[j]]] <- j
    }
    alloc.new[i] <- n_g
    groupings <- alloc.new
  }
  groupings
}

# Update permutation:
update_sigma <- function(sigma, k_rep, grouping, beta, delta, dist) {
  sigma_prop <- sigma
  ind <- sample(1:length(sigma), k_rep)
  sigma_prop[ind] <- sample(sigma_prop[ind])
  # accept/reject:
  a <- log_epa_prior(grouping, beta, delta, dist, sigma_prop) -
    log_epa_prior(grouping, beta, delta, dist, sigma)
  if (a > log(runif(1))){
    prop <- sigma_prop
  } else{
    prop <- sigma
  }
  prop
}

# Update mass parameter (beta):
update_beta <- function(beta, grouping, delta, dist, sigma) {
  beta_prop <- beta + rnorm(1, mean = 0, sd = .5)
  a <- log_epa_prior(grouping, beta_prop, delta, dist, sigma) + dgamma(beta_prop, 3, 5, log = T) -
    log_epa_prior(grouping, beta, delta, dist, sigma) - dgamma(beta, 3, 5, log = T)
  if (a > log(runif(1))) {
    prop <- beta_prop
  } else{
    prop <- beta
  }
  prop
}

# Update discount parameter (delta):
update_delta <- function(grouping, beta, delta, dist, sigma) {
  # Stepping in Slice Sampler:
  L <- 0
  U <- 1
  f_x <- exp(log_epa_prior(grouping, beta, delta, dist, sigma) + log_delta_prior(delta))
  y <- runif(1, 0, f_x)
  x <- runif(1, L, U)
  while (log(y) > log_epa_prior(grouping, beta, x, dist, sigma) + log_delta_prior(x)) {
    if (x < delta) {
      L = x
    } else{
      U = x
    }
    x <- runif(1, L, U)
  }
  x
}

# Function to run MCMC:(methods = aao or seq)
epa_mcmc <- function(N_i,
                     B = 10000,
                     thin = 1,
                     method = "aao",
                     prior.alpha,
                     g.a,
                     g.b,
                     dist) {
  # Fix n_i:
  non.na.ind <- !is.na(N_i[1, ])
  J <- ifelse(is.null(ncol(N_i)), length(N_i), ncol(N_i))
  n_i <- N_i[, !is.na(N_i[1, ])]
  if (is.null(nrow(n_i))) {
    n_i <- array(n_i, dim = c(1, length(n_i)))
  }
  K <- ifelse(is.null(nrow(n_i)), 1, nrow(n_i))
  # Set prior parameters:
  prior.alpha <- rep(prior.alpha, ncol(n_i))
  beta <- beta_sav <- .9
  delta <- delta_sav <- 0.01
  sigma <- 1:K
  k_rep <- floor(K / 2)

  ## Get initial values for alpha and gamma
  alpha <- array(NA, dim = c(1, ncol(n_i)))
  alpha[1, ] <- rep(1 / ncol(n_i), ncol(n_i))
  z <- alpha_to_z(log(alpha))
  z_1m <- sapply(z, log1mexp)
  gamma <- gamma_sav <- 5#rgamma(1, g.a, g.b)
  tables <- matrix(c(1, 1), nrow = 2)
  group_alloc <- 1
  groupings_sav <- array(NA, dim = c(B, nrow(n_i)))
  alpha_sav <- array(NA, dim = c(B, ncol(n_i)))
  theta_sav <- array(NA, dim = c(B, K, ncol(n_i)))
  # Possible partitions:
  p <- salso::enumerate.partitions(K)
  groupings <- p[sample(1:nrow(p), 1), ]
  # Run MCMC:
  # prior probs of partitions:

  # Run MCMC:
  for (b in 2:(B * thin)) {
    ## Update groupings:
    if (method == "aao") {
      groupings <- update_groupings_aao(p,
                                        beta,
                                        delta,
                                        dist,
                                        sigma,
                                        alpha,
                                        gamma,
                                        n_i,
                                        prior.alpha,
                                        g.a,
                                        g.b,
                                        K)
    } else{
      groupings <- update_groupings_seq(n_i,
                                        groupings,
                                        alpha,
                                        gamma,
                                        beta,
                                        delta,
                                        sigma,
                                        dist,
                                        prior.alpha,
                                        g.a,
                                        g.b,
                                        K)
    }

    ## Update permutation (sigma):
    sigma <- update_sigma(sigma, k_rep, groupings, beta, delta, dist)
    ## Update mass parameter:
    beta <- update_beta(beta, groupings, delta, dist, sigma)
    ## Update discount parameter (delta):
    delta <- update_delta(groupings, beta, delta, dist, sigma)
    # Update gamma/alpha:
    # Combine counts
    n_curr <- combine_counts(groupings, n_i)
    # Update gamma:
    gamma <- exp(gamma_update(z, z_1m, ncol(n_i), n_i = n_curr, nrow(n_curr) , log(gamma), g.a, g.b))

    # Update z:
    z <- update_z(z,
                  z_1m,
                  ncol(n_i),
                  n_i = n_curr,
                  nrow(n_curr),
                  log(gamma),
                  unique(prior.alpha))
    z_1m <- sapply(z, log1mexp)
    # Translate to alpha:
    alpha <- exp(alpha_map(z, z_1m))
    # Thin:
    if (b %% thin == 0) {
      alpha_sav[b / thin, ] <- alpha
      groupings_sav[b / thin, ] <- groupings
      gamma_sav[b / thin] <- gamma
      beta_sav[b / thin] <- beta
      delta_sav[b / thin] <- delta
      for (g in unique(groupings)) {
        ind <- which(groupings == g)
        theta_sav[b / thin, ind, ] <- matrix(
          rep(
            LaplacesDemon::rdirichlet(1, n_curr[g, ] + alpha * gamma),
            length(ind)
          ),
          nrow = length(ind),
          byrow = T
        )
      }
    }
  }

  alpha <- array(0, dim = c(B, J))
  alpha[, non.na.ind] <- alpha_sav
  theta <- array(0, dim = c(B, K, J))
  theta[, , non.na.ind] <- theta_sav
  return(
    list(
      groups = groupings_sav,
      alpha = alpha,
      gamma = gamma_sav,
      theta = theta,
      beta = beta_sav,
      delta = delta_sav
    )
  )
}

## Function for HAMC MCMC:
hamc_mcmc <- function(n_i, K, g.a, g.b, prior.alpha, B) {
  J <- sum(!is.na(n_i[1, ]))
  # Fix dimension of n_i if only one group:
  if (is.null(dim(n_i))) {
    n_i <- matrix(n_i, nrow = 1)
  }
  # initialize gamma and alpha (log scale)
  g <- c(log(1))
  alpha <- array(NA, dim = c(B + 500, J))
  theta <- array(0, dim = c(B + 500, K, length(n_i[1, ])))
  alpha[1, ] <- log(rep(1 / J, J))
  for (k in 1:K) {
    theta[1, k, !is.na(n_i[k, ])] <- LaplacesDemon::rdirichlet(1, exp(g[1] + alpha[1, ]) + n_i[k, ][!is.na(n_i[k, ])])
  }
  z <- z_1m <- psi <- array(NA, dim = c(B + 500, J - 1))
  z[1, ] <- alpha_to_z(alpha[1, ])
  psi[1, ] <- pbeta(exp(z[1, ]), .5, .5)
  z_1m[1, ] <- sapply(z[1, ], log1mexp)
  # Start MCMC:
  for (iter in 2:(B + 500)) {
    # Update Gamma:
    # Slice:
    g[iter] <- gamma_update(z[iter - 1, ], z_1m[iter - 1, ], J, n_i
                            , K, g[iter - 1], g.a, g.b)

    # Update alpha:
    ## Map to z(beta) variables:
    ## Update z values (does accept/reject with slice sampler):
    z[iter, ] <- update_z(z[iter - 1, ], z_1m[iter - 1, ], J, n_i = n_i
                          , K, g[iter], prior.alpha)
    psi[iter, ] <- pbeta(exp(z[iter, ]), .5, .5)
    z_1m[iter, ] <- sapply(z[iter, ], log1mexp)
    alpha[iter, ] <- alpha_map(z[iter, ], z_1m[iter, ])

    for (k in 1:K) {
      theta[iter, k, !is.na(n_i[k, ])] <- LaplacesDemon::rdirichlet(1, exp(g[iter] + alpha[iter, ]) + n_i[k, ][!is.na(n_i[k, ])])
    }
  }
  # Save row draws:
  alpha_draws <- array(0, dim = c(B, length(n_i[k, ])))
  alpha_draws[, !is.na(n_i[1, ])] <- exp(alpha[-c(1:500), ])
  gamma_draws <- exp(g[-c(1:500)])
  theta_draws <- theta[-c(1:500), , ]
  #print(i)

  return(list(
    alpha = alpha_draws,
    gamma = gamma_draws,
    theta = theta_draws
  ))
}
