## Log EPA prior:
log_epa_prior <- function(p, alpha, delta, dist, sigma){
  if(alpha < -delta | delta < 0 | delta >= 1){ # Set constraints
    -Inf
  }
  else {
    p_tm1 <- p[sigma[1]]
    q_tm1 <- 1
    count <- 2
    log_p_t <- log(1)
    for(i in sigma[2:length(p)]){
      if(p[i] %in% p_tm1){
        log_p_t[count] <- log(length(p_tm1) - delta * q_tm1)-log(alpha + length(p_tm1)) +
          log(sum(dist[sigma[which(p_tm1 == p[i])],i])) - log(sum(dist[sigma[1:length(p_tm1)],i]))
      } else{
        log_p_t[count] <- log(alpha + delta * q_tm1)-log(alpha + length(p_tm1))
      }
      count <- count + 1
      p_tm1 <- c(p_tm1, p[i])
      q_tm1 <- length(unique(p_tm1))
    }
    sum(log_p_t)
  }
}

log_delta_prior <- function(delta){
  log(.4*dbeta(delta, .05, 20) + .6*dbeta(delta, 2,2.5))
}
# Combine counts function:
combine_counts <- function(grouping, n_i){
  ids <- (unique(grouping))
  res <- lapply(ids, function(id) sort(c(1:length(grouping))[grouping == id]))
  names(res) <- 1:length(res)
  counts.list <- lapply(res, \(x) {
    if(length(x) == 1){
      n_i[x,]
    }
    else{colSums(n_i[x,])}})
  matrix(unlist(counts.list), ncol = ncol(n_i), byrow = T)
}
# Write full joint of gamma, alpha and counts
log_full_joint <- function(n, alpha, gamma, prior.alpha, g.a, g.b){
  a <- matrix(alpha, nrow = 1)
  pa <- matrix(prior.alpha, nrow = 1)
  ag <- alpha*gamma
  if(is.null(dim(n))){
    n <- matrix(n, nrow = 1)
  }
  n_ag <- t(apply(n,1,\(x) x + ag))
  LaplacesDemon::ddirichlet((a), (pa), log = T) + (g.a - 1)*log(gamma) - g.b*gamma +
    + nrow(n)*lgamma(gamma) + sum(-lgamma(rowSums(n_ag)) +
                                    rowSums(t(apply(lgamma(n_ag), 1, \(x){x - lgamma(ag)}))))
}
# Function to get posterior prob with given partition:
log_like_prob_group <- function(n_i, p, alpha, gamma, prior.alpha, g.a, g.b,
                                K){
  ids <- (unique(p))
  res <- lapply(ids, function(id) sort(c(1:K)[p == id]))
  names(res) <- 1:length(res)
  counts.list <- lapply(res, \(x) {
    if(length(x) == 1){
      n_i[x,]
    }
    else{colSums(n_i[x,])}})
  n_curr <- matrix(unlist(counts.list), ncol = ncol(n_i), byrow = T)
  log_full_joint(n_curr, alpha, gamma, prior.alpha, g.a, g.b)
}
# Function to update groupings using the all at once method:
update_groupings_aao <- function(p, beta, delta, dist, sigma, alpha, gamma, n_i, prior.alpha, g.a, g.b){
  log.prior.probs <- apply(p,1, \(x) log_epa_prior(x, beta, delta, dist, sigma))
  log.prior.probs <- log.prior.probs - max(log.prior.probs)
  ## Update groupings:
  log_probs <- apply(p, 1, \(x) log_like_prob_group(n_i, x,
                                                    alpha, gamma, prior.alpha, g.a, g.b, K))
  log_probs <- log_probs - max(log_probs)
  log_post_probs <- log_probs + log.prior.probs
  log_post_probs <- log_post_probs - max(log_post_probs)
  probs <- exp(log_post_probs)/sum(exp(log_post_probs))
  p[sample(1:52, 1, prob = probs),]
}

# Function for sequential group update:
update_groupings_seq <- function(n_i, groupings, alpha, gamma, beta, delta, sigma, dist, prior.alpha, g.a, g.b, K){
  count = 0
  for(i in sample(1:K)){
    alloc <- groupings[-i]
    # Combine n's:
    ids <- (unique(alloc))
    res <- lapply(ids, function(id) sort(c(1:K)[-i][alloc == id]))
    names(res) <- 1:length(res)
    counts.list <- lapply(res, \(x) {
      if(length(x) == 1){
        n_i[x,]
      }
      else{colSums(n_i[x,])}})
    n_curr <- matrix(unlist(counts.list), nrow = length(res), byrow = T)
    # Get log_post probabilities.
    log_probs <- log_lik <- log_prior <- c()
    for(j in 1:(length(ids) + 1)){
      n_use <- n_curr
      p_use <- numeric(K)
      p_use[i] <- j
      for(g in 1:length(ids)){
        p_use[res[[g]]] <- g
      }
      if(j != length(ids) + 1){
        n_use[j,] <- n_use[j,] + n_i[i,]
      } else{
        n_use <- rbind(n_use, n_i[i,])
      }
      log_lik[j] <- log_full_joint(n_use, alpha, gamma, prior.alpha, g.a, g.b)
      log_prior[j] <- log_epa_prior(p_use, beta, delta, dist, sigma)
    }
    log_lik <- log_lik - max(log_lik)
    log_prior <- log_prior - max(log_prior)
    log_probs <- log_lik + log_prior
    log_probs <- log_probs - max(log_probs)
    probs <- exp(log_probs)/sum(exp(log_probs))
    n_g <- sample(1:(length(ids) + 1), 1, prob = probs)
    alloc.new <- numeric(K)
    for(j in 1:length(ids)){
      alloc.new[res[[j]]] <- j
    }
    alloc.new[i] <- n_g
    groupings <- alloc.new
  }
  groupings
}

# Update permutation:
update_sigma <- function(sigma, k_rep, grouping, beta, delta, dist){
  sigma_prop <- sigma
  ind <- sample(1:length(sigma), k_rep)
  sigma_prop[ind] <- sample(sigma_prop[ind])
  # accept/reject:
  a <- log_epa_prior(grouping, beta, delta, dist, sigma_prop) -
    log_epa_prior(grouping, beta, delta, dist, sigma)
  if(a > log(runif(1))){
    prop <- sigma_prop
  } else{
    prop <- sigma
  }
  prop
}

# Update mass parameter (beta):
update_beta <- function(beta, grouping, delta, dist, sigma){
  beta_prop <- beta + rnorm(1, mean = 0, sd = .5)
  a <- log_epa_prior(grouping, beta_prop, delta, dist, sigma) + dgamma(beta_prop, 3, 5, log = T) -
    log_epa_prior(grouping, beta, delta, dist, sigma) - dgamma(beta, 3, 5, log = T)
  if(a > log(runif(1))){
    prop <- beta_prop
  } else{
    prop <- beta
  }
  prop
}

# Update discount parameter (delta):
update_delta <- function(grouping, beta, delta, dist, sigma){
  # Stepping in Slice Sampler:
  L <- 0
  U <- 1
  f_x <- exp(log_epa_prior(grouping, beta, delta, dist, sigma) + log_delta_prior(delta))
  y <- runif(1, 0, f_x)
  x <- runif(1, L, U)
  while(log(y) > log_epa_prior(grouping, beta, x, dist, sigma) + log_delta_prior(x)){
    if(x < delta){
      L = x
    } else{
      U = x
    }
    x <- runif(1, L, U)
  }
  x
}

# Function to run MCMC:(methods = aao or seq)
epa_mcmc <- function(N_i, B = 10000, thin = 1, method = "aao",
                     prior.alpha, g.a, g.b, dist){
  # Fix n_i:
  non.na.ind <- !is.na(N_i[1,])
  J <- ifelse(is.null(ncol(N_i)), length(N_i), ncol(N_i))
  n_i <- N_i[,!is.na(N_i[1,])]
  if(is.null(nrow(n_i))){n_i <- array(n_i, dim = c(1, length(n_i)))}
  K <- ifelse(is.null(nrow(n_i)), 1, nrow(n_i))
  # Set prior parameters:
  prior.alpha <- rep(prior.alpha, ncol(n_i))
  beta <- beta_sav <- .9
  delta <- delta_sav <- 0.01
  sigma <- 1:K
  k_rep <- floor(K/2)

  ## Get initial values for alpha and gamma
  alpha <- array(NA, dim = c(1, ncol(n_i)))
  alpha[1,] <- rep(1/ncol(n_i), ncol(n_i))#rdirichlet(1, prior.alpha)
  z <- alpha_to_z(log(alpha))
  z_1m <- sapply(z, log1mexp)
  gamma <- gamma_sav <- 5#rgamma(1, g.a, g.b)
  tables <- matrix(c(1,1), nrow = 2)
  group_alloc <- 1
  groupings_sav <- array(NA, dim = c(B, nrow(n_i)))
  alpha_sav <- array(NA, dim = c(B, ncol(n_i)))
  theta_sav <- array(NA, dim = c(B, K, ncol(n_i)))
  # Possible partitions:
  p <- salso::enumerate.partitions(K)
  groupings <- p[sample(1:nrow(p),1),]
  # Run MCMC:
  # prior probs of partitions:
  #log.prior.probs <- apply(p,1, \(x) log_epa_prior(x, beta[1], delta[1], dist, sigma[1,]))

  # Run MCMC:
  for(b in 2:(B*thin)){
    ## Update groupings:
    if(method == "aao"){
      groupings <- update_groupings_aao(p, beta, delta, dist, sigma, alpha, gamma, n_i,
                                        prior.alpha, g.a, g.b)
    }else{
      groupings <- update_groupings_seq(n_i, groupings, alpha, gamma,
                                        beta, delta, sigma, dist,
                                        prior.alpha, g.a, g.b, K)
    }

    ## Update permutation (sigma):
    sigma <- update_sigma(sigma, k_rep, groupings, beta,
                          delta, dist)
    ## Update mass parameter:
    beta <- update_beta(beta, groupings, delta, dist, sigma)
    ## Update discount parameter (delta):
    delta <- update_delta(groupings, beta, delta, dist, sigma)
    # Update gamma/alpha:
    # Combine counts
    n_curr <- combine_counts(groupings, n_i)
    # Update gamma:
    gamma <- exp(gamma_update(z, z_1m, ncol(n_i),  n_i = n_curr,
                              nrow(n_curr) , log(gamma), g.a, g.b))

    # Update z:
    z <- update_z(z, z_1m, ncol(n_i),  n_i = n_curr,
                  nrow(n_curr), log(gamma), unique(prior.alpha))
    z_1m <- sapply(z, log1mexp)
    # Translate to alpha:
    alpha <- exp(alpha_map(z, z_1m))
    # Thin:
    if(b %% thin == 0){
      alpha_sav[b/thin,] <- alpha
      groupings_sav[b/thin, ] <- groupings
      gamma_sav[b/thin] <- gamma
      beta_sav[b/thin] <- beta
      delta_sav[b/thin] <- delta
      for(g in unique(groupings)){
        ind <- which(groupings == g)
        theta_sav[b/thin,ind,] <- matrix(rep(LaplacesDemon::rdirichlet(1, n_curr[g,] + alpha*gamma), length(ind)),
                                         nrow = length(ind), byrow = T)
      }
    }
  }

  alpha <- array(0, dim = c(B, J))
  alpha[, non.na.ind] <- alpha_sav
  theta <- array(0, dim = c(B, K, J))
  theta[,,non.na.ind] <- theta_sav
  return(list(groups = groupings_sav, alpha = alpha, gamma = gamma_sav, beta = beta_sav, delta = delta_sav, theta = theta))
}
