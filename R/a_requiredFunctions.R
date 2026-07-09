# calculates log(1 - exp(x))
log1mexp <- function(x) {
  stopifnot(x <= 0)
  if (x > -0.693147) {
    # .693147 ~= log(2.0))
    out <- log(-expm1(x))
  } else if (x <= 0.0) {
    out <- log1p(-exp(x))
  }
  out
}
## Gets determinant for stick breaking transformation:

l.alpha.f.cond <- function(alpha, n_i, ga, prior.alpha) {
  k_obj <- c()
  if(!is.null(dim(n_i))){
    J <- ncol(n_i)
    K <- nrow(n_i)
  } else{
    J <- length(n_i)
    K <- 1
  }
  for (k in 1:K) {
    if(any(alpha == -Inf)){
      sum_obj <- lgamma(n_i[k, ][!is.na(n_i[k, ])] + exp(ga)* alpha) - lgamma(exp(ga)*
                                                                                alpha)
      sum_obj[which(alpha == -Inf)] <- 0
      k_obj[k] <- sum(sum_obj)
    } else{
      k_obj[k] <- sum(lgamma(n_i[k, ][!is.na(n_i[k, ])] + exp(ga)* alpha) - lgamma(exp(ga) *
                                                                                     alpha))
    }
  }

  ## Set up priors:
  log_prior <- sum((prior.alpha - 1) * log(alpha))
  res <- log_prior + sum(k_obj)
  return(res)
}

#Log Full-Conditional for gamma:
l.gamma.f.cond <- function(alpha, J, n_i, K, ga, g.a, g.b) {
  k_obj <- c()
  for (k in 1:K) {
    k_obj[k] <- lgamma(exp(ga)) + sum(lgamma(n_i[k, ][!is.na(n_i[k, ])] + exp(ga)*
                                               alpha)) - lgamma(sum(n_i[k, ], na.rm = T) + exp(ga)) - sum(lgamma(exp(ga)*
                                                                                                                   alpha))
  }
  ## gamma priors are g.a, g.b
  res <- (g.a - 1) * ga - exp(ga) * g.b + sum(k_obj)
  return(res)
}

## Functions for alpha slice sampler:
# Map z to alpha and alpha to z:
alpha_to_logits <- function(alpha) { # returns length one LESS than argument
  J <- length(alpha)
  log(alpha[1:(J-1)]) - log(alpha[J])
}

logits_to_alpha <- function(logits) { # returns length one MORE than argument
  elogits <- exp(logits)
  sumelogits <- sum(elogits)
  sump1 <- sumelogits + 1.0
  c(elogits, 1.0) / sump1
}

## marginal target for logits = logit(theta)
log_mtarg_logits <- function(logits, n_i, ga, prior.alpha) {
  alpha <- logits_to_alpha(logits)
  out0 <- l.alpha.f.cond(alpha, n_i, ga, prior.alpha)
  ljacob <- sum(log(alpha))
  out0 + ljacob
}

## Slice sampler functions:
slice_elliptical_mv_spherical <- function (z, sig, log_target) {
  nEvaluations <- 0
  k <- length(z)
  f <- function(z) {
    nEvaluations <<- nEvaluations + 1
    log_target(z)
  }
  fz <- f(z)$value
  stopifnot(fz > -Inf)
  y <- log(runif(1)) + fz
  nu <- rnorm(k, 0, sd = sig)
  twopi <- 2 * pi
  theta <- runif(1, 0, twopi)
  theta_min <- theta - twopi
  theta_max <- theta
  repeat {
    z1 <- z * cos(theta) + nu * sin(theta)
    fz1 <- f(z1)
    if (y < fz1$value) {
      return(list(z = z1, x = fz1$x, nEvaluations = nEvaluations))
    }
    if (theta < 0) {
      theta_min <- theta
    }
    else {
      theta_max <- theta
    }
    theta <- runif(1, theta_min, theta_max)
  }
}


slice_genelliptical_mv_logits <- function (x = NULL, # logits (transformed alpha)
                                           z = NULL, # transformed (standardized)
                                           log_target_logits, n_i, prior.alpha, gamma,
                                           mu, Sig, df, is_chol = FALSE,
                                           static_pseudo = FALSE # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
) {

  if (isTRUE(is_chol)) {
    SigL <- Sig
  }
  else {
    SigL <- t(chol(Sig))
  }

  if (isTRUE(static_pseudo)) { # was given z
    k <- length(z)
  } else { # was given x
    k <- length(x)
    z <- drop(forwardsolve(SigL, (x - mu)))
  }

  stopifnot(length(mu) == k)
  stopifnot(dim(Sig) == c(k, k))

  a <- 0.5 * (df + k)
  b <- 0.5 * (df + sum(z^2))
  s <- 1.0 / rgamma(1, shape = a, rate = b)

  lff <- function(zz) {

    xx <- drop(SigL %*% zz + mu)
    ssz <- sum(zz^2)
    out <- log_target_logits(xx, n_i, gamma, prior.alpha) + a * log1p(ssz/df)

    list(value = out, z = zz, x = xx)
  }

  out_ess <- slice_elliptical_mv_spherical(z = z, sig = sqrt(s), log_target = lff)

  out_ess
}


# Gamma slice sampler:
gamma_update <- function(alpha, J, n_i, k, ga, g.a, g.b) {
  ## CDF Transformation:
  y <-  l.gamma.f.cond(alpha, J, n_i, k, ga, g.a, g.b) + log(runif(1)) - log(pgamma(exp(ga), g.a, g.b))
  L <- 0
  R <- 1
  u <- runif(1, L, R)
  g <- log(qgamma(u, g.a, g.b))
  while (y >= l.gamma.f.cond(alpha, J, n_i, k, g, g.a, g.b) - log(pgamma(exp(g), g.a, g.b))) {
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

## Log EPA prior:
log_epa_prior <- function(p, alpha, delta, dist, sigma) {
  if (alpha < -delta | delta < 0 | delta >= 1) {
    # Set constraints
    -Inf
  } else {
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
  ag <- (alpha * gamma)
  if (is.null(dim(n))) {
    n <- matrix(n, nrow = 1)
  }
  n_ag <- t(apply(n, 1, \(x) x + ag))
  res <- LaplacesDemon::ddirichlet((a), (pa), log = T) + (g.a - 1) * log(gamma) - g.b *
    gamma +nrow(n) * lgamma((gamma)) + sum(-lgamma(rowSums(n_ag)) +
                                             rowSums(t(apply(
                                               lgamma(n_ag), 1, \(x) {
                                                 x - lgamma(ag)
                                               }
                                             ))))
  ifelse(is.nan(res), 0, res)
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
  p[sample(1:nrow(p), 1, prob = probs), ]
}

## Function to run the proportional likelihood of one group given a row of counts,
# partial_lik <- function(n, n_old, alpha,
#                             gamma){
#   res <- sum(lgamma(n + gamma*alpha)) - lgamma(sum(n) + gamma) +  lgamma(gamma) - sum(lgamma(gamma *alpha))
#   res_old <- sum(lgamma(n_old + gamma*alpha)) - lgamma(sum(n_old) + gamma) +  lgamma(gamma) - sum(lgamma(gamma *alpha))
#   #m <- max(res, res_old)
#   if(res < res_old){
#     -10e9
#   } else {
#     res + log1mexp(res_old - res)
#   }
# }
beta_log <- function(x){
  sum(lgamma(x))-lgamma(sum(x))
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
  for (i in 1:K) {
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
      p_use <- numeric(K)
      p_use[i] <- j
      for (g in 1:length(ids)) {
        p_use[res[[g]]] <- g
      }
      if (j != (length(ids) + 1)) {
        n_use <- n_curr[j, ] + n_i[i, ]
        log_lik[j] <- beta_log(n_use + gamma*alpha) - beta_log(n_curr[j,] + gamma*alpha) #+  lgamma(gamma) - sum(lgamma(gamma *alpha))
      } else{
        log_lik[j] <- beta_log(n_i[i,] + gamma*alpha) - beta_log(gamma*alpha)
      }

      log_prior[j] <- log_epa_prior(p_use, beta, delta, dist, sigma)
    }
    #log_lik <- log_full_joint(n_i[i,], alpha, gamma, prior.alpha, g.a, g.b)
    # log_lik <- log_lik - max(log_lik)
    # log_prior <- log_prior - max(log_prior)
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
  J <- sum(non.na.ind)#ifelse(is.null(ncol(N_i)), length(N_i), ncol(N_i))
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
  k_rep <- K#floor(K / 2)
  ## Get initial values for alpha and gamma
  alpha <- array(NA, dim = c(1, ncol(n_i)))
  alpha[1, ] <- rep(1 / ncol(n_i), ncol(n_i))
  x <- alpha_to_logits(alpha[1, ])

  gamma <- gamma_sav <- 5#rgamma(1, g.a, g.b)
  tables <- matrix(c(1, 1), nrow = 2)
  group_alloc <- 1
  groupings_sav <- array(NA, dim = c(B, nrow(n_i)))
  alpha_sav <- array(NA, dim = c(B, ncol(n_i)))
  theta_sav <- array(NA, dim = c(B, K, ncol(n_i)))
  # Possible partitions:
  if(K <= 13){
    p <- salso::enumerate.partitions(K)
    groupings <- p[sample(1:nrow(p), 1), ]
  } else {
    num_groups <- sample(1:K, 1)
    groupings <- sample(1:num_groups, K, replace = TRUE)
  }


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
    k <- length(unique(groupings))
    X <- n_curr
    Phat <- X / rowSums(X)[row(X)] # rowSums(X) is also n_vec
    phat_pool <- colSums(X) / sum(rowSums(X))
    phat_indep <- colMeans(Phat)
    nbar <- mean(rowSums(X))
    # Update Gamma:
    # Slice:
    gamma <- exp(gamma_update(alpha, J, n_curr
                            , k, log(gamma), g.a, g.b))

    ## Update mu/Sigma:
    ps_mix <- (gamma) * phat_pool + nbar * phat_indep + unique(prior.alpha)
    pseudo_scale <- 0.9
    pseudo_counts <- ps_mix * pseudo_scale # shape vector for Dirichlet pseudo-target

    ## Aitchison 1986, p 127 K-L approx for pseudo on logits
    mu <- digamma(pseudo_counts[-J]) - digamma(pseudo_counts[J])
    Sig <- matrix(trigamma(pseudo_counts[J]), nrow = J-1, ncol = J-1) + diag(trigamma(pseudo_counts[-J]))
    df <- k - 1
    # Update alpha:
    ## Map to z(beta) variables:
    ## Update z values (does accept/reject with slice sampler):
    tmp <- slice_genelliptical_mv_logits(x = x,
                                         z = NULL,
                                         log_target_logits = log_mtarg_logits, n_curr, prior.alpha, log(gamma),
                                         mu = mu, Sig = Sig, df = df, is_chol = FALSE,
                                         static_pseudo = FALSE # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
    )
    x <- tmp$x
    alpha <- logits_to_alpha(tmp$x)
    # Thin:
    if (b %% thin == 0) {
      alpha_sav[b / thin, ] <- alpha#[subset]
      groupings_sav[b / thin, ] <- groupings
      gamma_sav[b / thin] <- gamma
      beta_sav[b / thin] <- beta
      delta_sav[b / thin] <- delta
      for (g in unique(groupings)) {
        ind <- which(groupings == g)
        theta_sav[b / thin, ind, ] <- matrix(
          rep(
            LaplacesDemon::rdirichlet(1, n_curr[g, ] + alpha * gamma)[1,],#[1, subset],
            length(ind)
          ),
          nrow = length(ind),
          byrow = T
        )
      }
    }
  }
  J <- ifelse(is.null(ncol(N_i)), length(N_i), ncol(N_i))
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
hamc_mcmc <- function(N_i, K, g.a, g.b, prior.alpha, B) {
  J <- sum(!is.na(N_i[1, ]))

  # Fix dimension of n_i if only one group:
  if (is.null(dim(N_i))) {
    N_i <- matrix(N_i, nrow = 1)
  }
  # initialize gamma and alpha (log scale)
  g <- c(log(1))
  alpha <- array(NA, dim = c(B, J))
  theta <- array(0, dim = c(B, K, J))
  alpha[1, ] <- (rep(1 / J, J))
  for (k in 1:K) {
    theta[1, k, ] <- LaplacesDemon::rdirichlet(1, exp(g[1] + alpha[1, ]) + N_i[k, ][!is.na(N_i[k, ])])
  }
  x <- array(NA, dim = c(B, J - 1))
  x[1, ] <- alpha_to_logits(alpha[1, ])
  # Reorder by counts:
  n_i <- N_i[,!is.na(N_i[1, ])]
  X <- n_i
  Phat <- X / rowSums(X)[row(X)] # rowSums(X) is also n_vec
  phat_pool <- colSums(X) / sum(rowSums(X))
  phat_indep <- colMeans(Phat)
  nbar <- mean(rowSums(X))
  # Start MCMC:
  for (iter in 2:(B)) {
    # Update Gamma:
    # Slice:
    g[iter] <- gamma_update(alpha[iter -1,], J, n_i
                            , K, g[iter - 1], g.a, g.b)

    ## Update mu/Sigma:
    ps_mix <- exp(g[iter]) * phat_pool + nbar * phat_indep + unique(prior.alpha)
    pseudo_scale <- 0.9
    pseudo_counts <- ps_mix * pseudo_scale # shape vector for Dirichlet pseudo-target

    ## Aitchison 1986, p 127 K-L approx for pseudo on logits
    mu <- digamma(pseudo_counts[-J]) - digamma(pseudo_counts[J])
    Sig <- matrix(trigamma(pseudo_counts[J]), nrow = J-1, ncol = J-1) + diag(trigamma(pseudo_counts[-J]))
    df <- K - 1
    # Update alpha:
    ## Map to z(beta) variables:
    ## Update z values (does accept/reject with slice sampler):
    tmp <- slice_genelliptical_mv_logits(x = x[iter - 1,],
                                         z = NULL,
                                         log_target_logits = log_mtarg_logits, n_i, prior.alpha, g[iter],
                                         mu = mu, Sig = Sig, df = df, is_chol = FALSE,
                                         static_pseudo = FALSE # could mu and Sig change from iteration to iteration? If not, we can build the chain on Z and save computation
    )
    x[iter,] <- tmp$x
    alpha[iter,] <- logits_to_alpha(tmp$x)
    for (k in 1:K) {
      theta[iter, k, ] <- LaplacesDemon::rdirichlet(1, exp(g[iter])*alpha[iter, ] + n_i[k, ])
    }
  }
  # Save row draws:
  alpha_draws <- array(0, dim = c(B, length(N_i[1, ])))
  theta_draws <- array(0, dim = c(B, K, length(N_i[1, ])))
  alpha_draws[, !is.na(N_i[1, ])] <- (alpha)
  gamma_draws <- exp(g)
  theta_draws[ , , !is.na(N_i[1, ])] <- theta[, , ]
  #print(i)

  return(list(
    alpha = alpha_draws,
    gamma = gamma_draws,
    theta = theta_draws
  ))
}
