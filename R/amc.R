#' AMC model
#'
#' @description
#' This function fits an indpendent AMC model on each row of the transition matrix and outputs posterior draws of the transition probabilities.
#'
#'
#' @param N A two or three dimensional array of counts (KxIxJ) or (IxJ) where K is the number of subgroups, I is the number of transitory states, J is the total number of states.
#' @param B Number of posterior draws to return (default is 10000)
#' @param xsi the a vector of length I where xsi is the prior value for the uniform Dirichlet prior for the ith row. The higher xsi will lead to stronger priors towards uniform transition probabilities.
#'
#' @return A three dimensional array that is BxIxJ of the posterior draws of the transition matrices.
#' @export
#' @examples
#'


amc <- function(N, B = 10000, xsi = NULL) {
  # Check compatability of N:
  d <- dim(N)
  if (!is.array(N)) {
    stop(paste("N must be an array."))
  }
  if (length(d) < 2) {
    stop(paste("N must be a two orthree dimensional array."))
  }
  if (length(d) == 3) {
    N_new <- array(NA, dim = c(d[2], d[3]))
    for (i in 1:d[2]) {
      N_new[i, ] <- colSums(N[, i, ], na.rm = TRUE)
    }
    N <- N_new
  }

  I <- nrow(N)
  J <- ncol(N)

  if (!is.null(xsi)) {
    if (xsi <= 0) {
      stop(paste("xsi must be positive"))
    }
    if (length(xsi) == 1) {
      xsi <- rep(xsi, I)
    }
  }
  if (is.null(xsi)) {
    xsi <- c()
    for (i in 1:I) {
      xsi[i] <- 1 / sum(!is.na(N[i, ]))
    }
  }

  ## Get posterior draws:
  theta <- array(0, dim = c(B, I, J))
  for (i in 1:I) {
    theta[, i, !is.na(N[i, ])] <- LaplacesDemon::rdirichlet(B, N[i, !is.na(N[i, ])] + rep(xsi[i], sum(!is.na(N[i, ]))))
  }

  return(theta)
}
