#' HAMC model
#'
#' @description
#' This function fits an indpendent HAMC model on each row of the transition matrix and outputs posterior draws of the two hyper parameters (alpha and gamma) and the transition probabilities.
#'
#'
#' @param N A two or three dimensional array (KxIxJ) where K is the number of subgroups, I is the number of transitory states, J is the total number of states.
#' @param B Number of posterior draws to get (default is 10000)
#' @param xsi the a vector of length I where xsi is the prior value for the uniform dirichlet prior.
#' @param g.a The prior shape parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.
#' @param g.b The prior rate parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.
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
