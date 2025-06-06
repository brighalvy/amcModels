% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epa_hamc.R
\name{epa_hamc}
\alias{epa_hamc}
\title{HAMC-EPA model}
\usage{
epa_hamc(
  N,
  method = "aao",
  B = 10000,
  thin = 1,
  prior.alpha = NULL,
  g.a = NULL,
  g.b = NULL,
  dist = NULL,
  MCMC.cores = 1,
  burnin = 500
)
}
\arguments{
\item{N}{A three dimensional array of counts (KxIxJ) where K is the number of subgroups, I is the number of transitory states, J is the total number of states. Include NA's in matrix for impossible transitions.}

\item{method}{The MCMC approach to obtaining posterior draws of the random clusters. Options "aao" or "seq". "aao" is an "all-at-once" approach in that it evaluates the posterior probability of each possible partition. "seq" does a sequential allocation of clustering given the previous clustering draw. When there are a small number of subgroups "aao" may be prefered however as the number of subgroups increases the computation time increases dramatically.}

\item{B}{Number of desired posterior estimates}

\item{thin}{Thinning of MCMC draws, i.e. thin = 5, means every 5th draw will be returned.}

\item{prior.alpha}{The value to be included in the prior for the Dirichlet on alpha, will be uniform. Default is 1/(# of possible transition for each row).}

\item{g.a}{The prior shape parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.}

\item{g.b}{The prior rate parameter of the gamma distribution for the gamma parameter. Default uses a mean of 20 and standard deviation of 18.}

\item{dist}{A KxK symmetric matrix that includes the prior information on the distance between subgroups. Defaults to each group is sequentially one from the next in line.}

\item{MCMC.cores}{The number of cores desired to use to run, if greater than 1 parallel processing will be used. Default is 1 (no parallelization).}

\item{burnin}{The number of burn in iterations to be removed from the posterior draws. Default is 500.}
}
\value{
A list with the following:
\itemize{
\item groups - The posterior draws of the clustering from each MCMC draw. (BxIxK)
\item alpha - Posterior draws of the alpha parameter. (BxIxJ)
\item gamma - Posterior draws of the gamma parameter. (BxI)
\item theta - Posterior draws of the transition probabilities. (BxKxIxJ)
}
}
\description{
This function fits an indpendent HAMC-EPA model on each row of the transition matrix and outputs posterior draws of clustering, alpha, gamma, and the transition probabilities.
}
