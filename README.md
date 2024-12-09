# Absorbing Markov Chain Models (amcModels)

The purpose of this package is to implement the methodology in "Hierarchical Absorbing Markov Chain Models for Pooling Sparse Combat Sequence Data". This paper is currently in the editing stages for review. A reference will be added once it becomes available.

## Intended Functionality:

The main functionality of this package is to fit one of the three absorbing Markov chain models. One function `amc` will require a matrix of counts and a desired vector for the Dirichlet prior and will return posterior draws of the transition matrix, it will require the `LaplacesDemon` package. The second function will fit the model with a hierarchical structure. This method requires an array of counts split by subgroup, a desired vector for the Dirichlet prior on $\boldsymbol{\alpha}_i$ and the prior parameters for $\gamma_i$. This function will be names `hamc`. It will return posterior draws of the transition matrices, $\boldsymbol{\alpha}_i$ and $\gamma_i$. The third function will be the extension using the Ewen's Pitman attraction. This function will be called `epa_hamc`. It requires all of the same arguments as `hamc` in addition it now requires a distance matrix that will be a $K \times K$ matrix outlining your prior belief of similarities between subgroups. It will output the same as `hamc` but will add in addition the posterior draws of the potential clusterings of the subgroups. In addition a function `absmarkovchain` will be added to be able to simulate data from a given transition matrix.

## Installation Instructions:

-   If `devtools` is not installed, install devtools.

-   Run the following line of code: `devtools::install_github("brighalvy/amcModels")`

## Small Example:

We will simulate sequence data from 2 groups with 2 transitory states and 1 absorbing state. After simulating the data we will fit each model on the created data.

```{r}
library(amcModels)

# Set transition matrices:
P1 <- matrix(c(1/3, 1/3, 1/3,
               .3, .5, .2,
               .5, .4, .1), nrow = 3, ncol = 3, byrow = TRUE)
P2 <- matrix(c(.3, .4, .4,
               .5, .3, .2,
               .45, .1, .45))
               
# Simulate the data
n <- 100

## Simulate data starting in state 1:
data_1_1 <- absmarkovchain(1, P1, 100)
data_2_1 <- absmarkovchain(1, P2, 100)

## Simulate data starting in state 2:
data_1_2 <- absmarkovchain(2, P1, 100)
data_2_2 <- absmarkovchain(2, P2, 100)

# Create count array
N <- array(0, dim = c(6, 3, 5))

## Function to apply to each sequence:
make_counts <- function(seq, N_group){
  for(x in 2:length(seq)){
    row.ref <- seq[x - 1]
    col.ref <- seq[x]
    N_group[row.ref, col.ref] <- N_group[row.ref, col.ref] + 1
  }
  return(N_group)
}

for(i in 1:100){
  ## Group 1:
  N[1, , ] <- make_counts(data_1_1[[i]], N[1, , ])
  N[1, , ] <- make_counts(data_1_2[[i]], N[1, , ])
  
  ## Group 2:
  N[2, , ] <- make_counts(data_2_1[[i]], N[2, , ])
  N[2, , ] <- make_counts(data_2_2[[i]], N[2, , ])

}

# Fit amc model:
fit_amc <- amc(N, xsi = 1/3)

# Fit hamc model:
fit_hamc <- hamc(N, MCMC.cores = 1)

# Fit epa-hamc model:
fit_epa <- epa_hamc(N, method = "seq", thin = 5, MCMC.cores = 1)
```
