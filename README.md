\section{Description}

One approach to analyzing the progression of discrete events through time is using Markov chain models. Sometimes these sequences can have what is called absorbing states which are states that end a sequence of events, sometimes known as terminal states thus we can use an absorbing Markov chain. A Markov chain is equipped with a transition matrix that contains transition probabilities between states. Our goal is to model this transition matrix. There are three different approach we will take. The first is using the dataset as a whole and fitting an independent Dirichlet prior on each row of the transition matrix. In certain scenarios however there exist subgroups ($k$) that effect the transition probabilities and it would be useful to analyze the differences in the transition probabilities. To do so we can add a hierarchical structure as follows:

\begin{align}
    \{x_{gl}\}_{g = k} \mid \mathbf{P}^{(k)} &\sim \text{Absorbing Markov Chain}(\mathbf{P}^{(k)}), \text{ where } k \in \{1, 2, ..., K\},\\
    \mathbf{p}_i^{(k)} \mid \boldsymbol{\alpha}_i, \gamma_i &\sim \text{Dirichlet}(\boldsymbol{\alpha}_i\gamma_i), \\
    \boldsymbol{\alpha}_i &\sim \text{Dirichlet}(\boldsymbol{\xi}), \gamma_i \sim \text{Gamma}(\tau, \eta).
\end{align}

The two methods, with and without hierarchy, work generally well in most cases however in some instances the data may be spread sparsely over the subgroups. To combat this you can add random clustering to the hierarchical model to then allow for some of the subgroups to be combined in their inference. To do so we will use what is called the Ewen's Pitman Attraction distribution to provide groupings and then use the hierarchical structure on the new groupings..

\section{Functionality}

The main functionality of this package is to fit one of the three models described above. One function \verb|amc| will require a matrix of counts and a desired vector for the Dirichlet prior and will return posterior draws of the transition matrix, it will require the \verb|LaplacesDemon| package.

The second function will fit the model with a hierarchical structure. This method requires an array of counts split by subgroup, a desired vector for the Dirichlet prior on $\boldsymbol{\alpha}_i$ and the prior parameters for $\gamma_i$. This function will be names \verb|hamc|. It will return posterior draws of the transition matrices, $\boldsymbol{\alpha}_i$ and $\gamma_i$.

The third function will be the extension using the Ewen's Pitman attraction. This function will be called \verb|epa_hamc|. It requires all of the same arguments as \verb|hamc| in addition it now requires a distance matrix that will be a $K \times K$ matrix outlining your prior belief of similarities between subgroups. It will output the same as \verb|hamc| but will add in addition the posterior draws of the potential clusterings of the subgroups.
