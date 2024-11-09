# Absorbing Markov Chain Models (amcModels)

## Intended Functionality:

The main functionality of this package is to fit one of the three models described above. One function \verb|amc| will require a matrix of counts and a desired vector for the Dirichlet prior and will return posterior draws of the transition matrix, it will require the \verb|LaplacesDemon| package.

The second function will fit the model with a hierarchical structure. This method requires an array of counts split by subgroup, a desired vector for the Dirichlet prior on $\boldsymbol{\alpha}_i$ and the prior parameters for $\gamma_i$. This function will be names \verb|hamc|. It will return posterior draws of the transition matrices, $\boldsymbol{\alpha}_i$ and $\gamma_i$.

The third function will be the extension using the Ewen's Pitman attraction. This function will be called \verb|epa_hamc|. It requires all of the same arguments as \verb|hamc| in addition it now requires a distance matrix that will be a $K \times K$ matrix outlining your prior belief of similarities between subgroups. It will output the same as \verb|hamc| but will add in addition the posterior draws of the potential clusterings of the subgroups.


## What is left?

I still need to create the amc and hamc functions, these will be built on the required functions already created. Further I will need to optimize my function implementations to speed up the MCMC process.
