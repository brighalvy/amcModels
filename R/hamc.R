#' HAMC model
#'
#' @description
#' This function fits an indpendent HAMC model on each row of the transition matrix and outputs posterior draws of the two hyper parameters (alpha and gamma) and the transition probabilities.
#'
#'
#' @param N A two or three dimensional array (KxIxJ) where K is the number of subgroups, I is the number of transitory states, J is the total number of states.
#' @param B Number of posterior draws to get (default is 10000)
#' @param prior.alpha the a vector of length I where xsi is the prior value for the uniform dirichlet prior.
#' @param g.a The prior shape parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.
#' @param g.b The prior rate parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.
#' @param MCMC.cores the number of cores to use for parralel processing, the default is 1.
#'
#' @return A list with the following:
#'  \itemize{
#'    \item alpha - Posterior draws of the alpha parameter. (BxIxJ)
#'    \item gamma - Posterior draws of the gamma parameter. (BxI)
#'    \item theta - Posterior draws of the transition probabilities. (BxKxIxJ)
#'  }
#' @export
#' @examples
#'


hamc <- function(N, B = 10000, prior.alpha = NULL, g.a = NULL, g.b = NULL, MCMC.cores = 1){
  # Check compatability of N:
  if (!is.array(N)) {
    stop(paste("N must be an array."))
  }
  if (length(dim(N)) < 3) {
    stop(paste("N must be a three dimensional array."))
  }
  # Check B, thin, prior.alpha, g.a, g.b are numeric if input is supplied:
  if (!is.numeric(B)) {
    stop(paste("B must be numeric"))
  }
  if (!is.null(g.a) & !is.numeric(g.a)) {
    stop(paste("g.a must be numeric"))
  }
  if (!is.null(g.b) & !is.numeric(g.b)) {
    stop(paste("g.b must be numeric"))
  }
  # Make sure prior parameters are positive:
  if (!is.null(prior.alpha)) {
    if(prior.alpha <= 0){
      stop(paste("prior.alpha must be positive"))
    }
  }
  if(is.null(prior.alpha)){
    prior.alpha <- c()
    for(i in 1:dim(N)[2]){
      prior.alpha[i] <- 1/sum(!is.na(N[1, i, ]))
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
  if(length(unique(prior.alpha)) == 1){
    prior.alpha <- unique(prior.alpha)
    fit <- parallel::mclapply(count.list,
                              \(x) hamc_mcmc(x, K, g.a, g.b, prior.alpha, B),
                              mc.cores = MCMC.cores) # n_i, K, g.a, g.b, prior.alpha, B
  } else{
    for(i in 1:dim(N)[2]){
      count.list[[i]] <- cbind(count.list[[i]], prior.alpha[i])
    }
    fit <- parallel::mclapply(count.list,
                              \(x) hamc_mcmc(x[, 1:(ncol(x) - 1)], K, g.a, g.b, x[1, ncol(x)], B),
                              mc.cores = MCMC.cores)
  }


  ## Create output list options:
  I <- dim(N)[2]
  J <- dim(N)[3]
  alpha <- array(0, dim = c(B, I, J))
  gamma <- array(0, dim = c(B, I))
  theta <- array(0, dim = c(B, K, I, J))
  for (i in 1:I) {
    alpha[, i, ] <- fit[[i]][[1]]
    gamma[, i] <- fit[[i]][[2]]
    theta[, , i, ] <- fit[[i]][[3]]
  }
  output <- list(
    alpha = alpha,
    gamma = gamma,
    theta = theta
  )


  return(output)
}
