#' HAMC-EPA model
#'
#' @description
#' This function fits an indpendent HAMC-EPA model on each row of the transition matrix and outputs posterior draws of clustering, alpha, gamma, and the transition probabilities.
#'
#'
#' @param N A three dimensional array (KxIxJ) where K is the number of subgroups, I is the number of transitory states, J is the total number of states.
#' @param method The MCMC approach to obtaining posterior draws of the random clusters. Options "aao" or "seq". "aao" is an "all-at-once" approach in that it evaluates the posterior probability of each possible partition. "seq" does a sequential allocation of clustering given the previous clustering draw. When there are a small number of subgroups "aao" may be prefered however as the number of subgroups increases the computation time increases dramatically.
#' @param B Number of desired posterior estimates
#' @param thin Thinning of MCMC draws, i.e. thin = 5, means every 5th draw will be returned.
#' @param prior.alpha The value to be included in the prior for the Dirichlet on alpha, will be uniform. Default is 1/(# of possible transition for each row).
#' @param g.a The prior shape parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.
#' @param g.b The prior rate parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.
#' @param dist A KxK symmetric matrix that includes the prior information on the distance between subgroups. Defaults to each group is sequentially one from the next in line.
#' @param MCMC.cores The number of cores desired to use to run, if greater than 1 parallel processing will be used. Default is 1 (no parallelization).
#'
#' @return A list with the following:
#'  \itemize{
#'    \item groups - The posterior draws of the clustering from each MCMC draw. (BxIxK)
#'    \item alpha - Posterior draws of the alpha parameter. (BxIxJ)
#'    \item gamma - Posterior draws of the gamma parameter. (BxI)
#'    \item theta - Posterior draws of the transition probabilities. (BxKxIxJ)
#'  }
#' @export
#' @examples
#'

# Write function:
epa_hamc <- function(N,
                     method = "aao",
                     B = 10000,
                     thin = 1,
                     prior.alpha = NULL,
                     g.a = NULL,
                     g.b = NULL,
                     dist = NULL,
                     MCMC.cores = 1) {
  # Check compatability of N:
  if (!is.array(N)) {
    stop(paste("N must be an array."))
  }
  if (length(dim(N)) < 3) {
    stop(paste("N must be a three dimensional array."))
  }
  # Check method argument:
  if (!(method %in% c("aao", "seq"))) {
    stop(paste("Invalid method argument"))
  }
  # Check B, thin, prior.alpha, g.a, g.b are numeric if input is supplied:
  if (!is.numeric(B)) {
    stop(paste("B must be numeric"))
  }
  if (!is.numeric(thin)) {
    stop(paste("thin must be numeric"))
  }
  if (!is.null(g.a) & !is.numeric(g.a)) {
    stop(paste("g.a must be numeric"))
  }
  if (!is.null(g.b) & !is.numeric(g.b)) {
    stop(paste("g.b must be numeric"))
  }
  # Make sure prior parameters are positive:
  if (!is.null(prior.alpha)) {
    if (prior.alpha <= 0) {
      stop(paste("prior.alpha must be positive"))
    }
  }
  if (is.null(prior.alpha)) {
    prior.alpha <- c()
    for (i in 1:dim(N)[2]) {
      prior.alpha[i] <- 1 / sum(!is.na(N[1, i, ]))
    }
  }
  if (!is.null(g.a)) {
    if (g.a <= 0) {
      stop(paste("g.a must be positive"))
    }
  }
  if (!is.null(g.b)) {
    if (g.b <= 0) {
      stop(paste("g.b must be positive"))
    }
  }
  # Check distance matrix and correct dimensions:
  K <- dim(N)[1]
  if (!is.null(dist)) {
    if (!is.matrix(dist)) {
      stop("dist must be a matrix.")
    }
    if (dim(dist)[1] != K) {
      stop("Wrong dimensions for distance matrix.")
    }
  }
  # Check dimensions
  if (is.null(dist)) {
    dist <- matrix(1, nrow = K, ncol = K)
    for (k in 1:K) {
      dist[c(1:K)[-k], k] <- abs(c(1:K)[-k] - k)
    }
    dist[K, K] <- 1
    # Convert to similarity matrix for EPA model:
    dist <- dist ^ (-.5)
  } else {
    # Convert to similarity matrix for EPA model:
    dist <- dist ^ (-.5)
  }
  ## Populate g.a and g.b if missing:
  if (is.null(g.a) | is.null(g.b)) {
    mu <- 20
    sd <- 18
    g.b <- mu / sd ^ 2
    g.a <- mu * g.b
  }

  # Make the data into a list:
  count.list <- list()
  for (i in 1:dim(N)[2]) {
    count.list[[i]] <- N[, i, ]
  }

  ## Fit model:
  ## Check configuration of prior.alpha if it varies from row to row:
  if (length(unique(prior.alpha)) == 1) {
    prior.alpha <- unique(prior.alpha)
    fit <- parallel::mclapply(count.list,
                              \(x) epa_mcmc(x, B, thin, method, prior.alpha, g.a, g.b, dist),
                              mc.cores = MCMC.cores)
  } else{
    for (i in 1:dim(N)[2]) {
      count.list[[i]] <- cbind(count.list[[i]], prior.alpha[i])
    }
    fit <- parallel::mclapply(count.list,
                              \(x) epa_mcmc(x[, 1:(ncol(x) - 1)], B, thin, method, x[1, ncol(x)], g.a, g.b, dist),
                              mc.cores = MCMC.cores)
  }


  ## Create output list options:
  I <- dim(N)[2]
  J <- dim(N)[3]
  groupings <- array(NA, dim = c(B, I, K))
  alpha <- array(0, dim = c(B, I, J))
  gamma <- array(0, dim = c(B, I))
  theta <- array(0, dim = c(B, K, I, J))
  for (i in 1:I) {
    groupings[, i, ] <- fit[[i]][[1]]
    alpha[, i, ] <- fit[[i]][[2]]
    gamma[, i] <- fit[[i]][[3]]
    theta[, , i, ] <- fit[[i]][[4]]
  }
  output <- list(
    groups = groupings,
    alpha = alpha,
    gamma = gamma,
    theta = theta
  )


  return(output)
}
