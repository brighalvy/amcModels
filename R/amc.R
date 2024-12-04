#' AMC model
#'
#' @description
#' This function fits an indpendent AMC model on each row of the transition matrix and outputs posterior draws of the transition probabilities.
#'
#'
#' @param N A two or three dimensional array (KxIxJ) where K is the number of subgroups, I is the number of transitory states, J is the total number of states.
#' @param B Number of posterior draws to get (default is 10000)
#' @param xsi the value to be included as the prior parameter in the uniform dirichlet prior for each row.
#'
#' @return A three dimensional array
#' @export
#' @examples
#'


amc <- function(N, B = 10000, xsi){
  # Check compatability of N:
  if (!is.array(N)) {
    stop(paste("N must be an array."))
  }
  if (length(dim(N)) < 3) {
    stop(paste("N must be a three dimensional array."))
  }
}
