#' Absorbing Markov Chain
#'
#' @description
#' This function returns n sequences of events from an absorbing Markov chain.
#'
#' @param x0 initial state, integer valued.
#' @param P transition probability matrix (IxJ) where the absorbing states are the last J-I rows. I = number of transitory states, J = total number of states.
#' @param n the number of sequences.
#'
#' @return A list where each element is a sequence of events from P.
#' @export
#' @examples
#'


absmarkovchain <- function(x0, P, n) {
  # Dimensions
  I <- nrow(P)
  J <- ncol(P)

  ## Set up run:
  sequences <- list()
  for (j in 1:n) {
    abs_states <- (I + 1):J
    states <- 1:J
    trans <- 1:I
    chain <- x0
    i = 1
    while ((chain[i] %in% trans)) {
      i = i + 1
      chain[i] <- sample(states, size = 1, prob = P[chain[i - 1], ])
    }
    sequences[[j]] <- chain
  }

  return(sequences)
}
