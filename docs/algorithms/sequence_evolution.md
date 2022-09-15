
# Sequence evolution models

Sequence evolution models describe the stochastic ways sequences change over time do to mutations.
If we consider just a single site and denote the probability of finding that site in state $i$ as $p_i$, these models take the form
$$\frac{d p_i}{dt} = \sum_j Q_{ij} p_j$$
where the matrix $Q_{ij}$ describes the rate at which the mutation changes the site from state $j$ to $i$.
This is a linear equation and it is solved by
$$p_i(t) = \left( e^{\mathbf{Q} t} \right)\_{i,j} p_j(0)$$
where $\left(e^{\mathbf{Q} t}\right)_{ij}$ is the $i,j$ element of the matrix exponential of $\mathbf{Q} t$.
This matrix exponential is central to sequence evolution and one of the most computationally demanding parts.

Since the matrix $Q_{ij}$ describes the evolution of a probability distribution, it needs to conserve probability.
In other words the sum $\sum_i p_i=1$ and
$$\sum_i \frac{d p_i}{dt} = 0 = \sum_{ij} Q_{ij} p_j$$
Since this has to be valid for any $p_j$ (take for example the case where $p_j=1$ and $p_i =0$ for $i\neq j$), we need $\sum_{i} Q_{ij} = 0$ which is typically enforced by requiring $Q_{jj} = - \sum_{i\neq j} Q_{ij}$.

This has the additional consequence that there has to be at least one vanishing eigenvalue.
Typically, this is exactly one eigenvalue but there could be special non-reducible models that have more than one (we have never considered such cases).
We adopt the convention that the vanishing eigenvalue is $\lambda_n$ where $n$ is the size of the alphabet.
The corresponding right eigenvalue $v^n_i$ is the equilibrium distribution $\pi_i$.
$$0 = \frac{d \pi_i}{dt} = \sum_j Q_{ij} \pi_j$$
Note that the corresponding left eigenvalue is a vector $\bar{v}^n_i$ with just $\bar{v}^n_i = 1$ at every entry.
$$0 = \sum_i \bar{v}^n_i Q_{ij} = \sum_i  Q_{ij} = 0$$
The left and right eigenvectors are ortho-normal, we can further deduce that
$$\sum_i \bar{v}^{n}_i v^m_i = \sum_i v^m_i = 0 \quad \mathrm{m\neq n}$$
that all right eigenvectors with non-zero eigenvalues have to sum to zero.

## Average rate
Traditionally, sequence evolution models are normalized such that the average rate equals one, that is on a branch of length one you observe on average one mutation.
The rate of mutating away from site $i$ is $Q_{ii} = -\sum_j Q_{ji}$ and the average rate is thus $-\sum_{i} Q_{ii}\pi_i = 1$.
If the average rate calculated such differs from 1, $Q_{ij}$ can be rescaled accordingly (we can't rescale $\pi_i$ as they sum to one).



## General time reversible models

Most models that are used in phylogenetics are so called time reversible models with the property that in equilibrium there are no net fluxes, that is that the flux from state $i$ to $j$ equals the flux from $j$ to $i$:
$$Q_{ij} \pi_j = Q_{ji}\pi_i$$
A simple and general way to enforce this is to parameterize the rate matrix as
$$Q_{ij} = \pi_i W_{ij} \quad \mathrm{i\neq j}$$
where $W_{ij}$ is a symmetric matrix.
The diagonal $Q_{jj} = -\sum_i Q_{ij}$ as before.

### Eigenvector calculations for general GTR models

Calculating left and right eigenvectors of a general matrix is numerically challenging.
But for GTR models, this problem can be reduced to the calculation of eigenvectors of a symmetric matrix.
Instead of calculating eigenvectors of $Q_{ij}$, we define the diagonal matrix $D_{ij} = \delta_{ij}\sqrt{\pi_j}$ and consider
$$\tilde{\mathbf{Q}} = \mathbf{D^{-1} Q D} = \sqrt{\pi_i} W_{ij}\sqrt{\pi_j}$$
We can now calculate eigenvalues and eigenvectors $\tilde{\mathbf{Q}}$, for which being a symmetric matrix the right and left eigenvectors agree and standard routines such as $\texttt{eigh}$ can be used.
An eigenvector $w^k_i$ to $\tilde{\mathbf{Q}}$ obeys
$$\mathbf{D} \mathbf{\tilde{Q}}  \mathbf{D}^{-1} \mathbf{D} w^k_i = \mathbf{D} \lambda_k  w^k_i = \mathbf{Q}  \mathbf{D} w^k_i$$
and $\mathbf{D} w^k_i$ is thus a right eigenvector to $\mathbf{Q}$ with eigenvector $\lambda_k$.
Similarly, the left eigenvectors are $\mathbf{D}^{-1} w^k_i$.


## Special models

### Jukes Cantor (JC69)

The Jukes-Cantor model is the simplest possible substitution model with an equal rate between any pair of states $Q_{ij} = 1/(n-1)$ where $n$ is the alphabet size.
In this case, $\lambda_i = (n+1)/n$ for $i=1\ldots, n-1$ and $\lambda_n=0$.
The eigenvector $v^n_j=1$ and the remaining $n-1$ are degenerate and span the orthogonal space.
Left and right eigenvectors agree since the model is symmetric.

The sequence evolution of this model is particularly simple
$$P(s|t) = (1-e^{-(n+1)t/n})/n + e^{-(n+1)t/n}P(s|0)/n$$
The advantage of this model is that we only need to calculate $n$ numbers rather than calculate multiple matrix-vector dot-products.


### HKY model

The HKY model is a fairly general model that allows for uneven base frequencies and transition/transversion biases but can still be calculated analytically.
However, evolving the state vector still involves a matrix product.
It is thus not clear whether this is worth implementing, the only aspect saved would be the one-time calculation of eigenvectors.



### Inference of GTR models

The most likely GTR models can be inferred efficiently from counts of observed substitutions.
This method was described in [Puller et al](https://doi.org/10.1093/ve/veaa066).

### Site specific models

Different sites in the sequence might differ both in the rate $\mu$ at which they change and in the relative rate of the changes.
The former can be implemented by storing a vector $\mu^a$ for each site $a$ instead of just one constant rate.
The eigenvalues and eigenvectors would be shared among all sites.

If instead the equilibrium frequencies $\pi$ vary from site to site, then eigenvalues and eigenvectors change along the sequence.
This implies that the matrix $e^{Q^at}$ is site specific which increases the computational burden.

## Organization of sequence data and models

Traditionally, sequence data comes in the form of an alignment with one sequence per individually.
At each node, the sequence would then be stored as a string of characters or a matrix describing the likelihood of each sequence state
```
seq A   N   C   G   ...
pos 0   1   2   3   ...
A   1   1   0   0   ...
C   0   1   1   0
G   0   1   0   1
T   0   1   0   0
```
These profiles consist of the column vectors that enter the sequence evolution models.
Storing these profiles can be quite costly, as they require `nL` floating point numbers.
For a large data set (say 10k SARS-CoV-2 sequences) this would require $10k\dot 30k\dot 5 \dot 8b=600Mb$.
In fact, this would need to be further multiplied by 2 to account for internal nodes.
Current TreeTime also keeps a few intermediate results around, further increasing the memory footprint.
Memory should therefore be a serious consideration.

One common way in which sequences can be compressed is by lumping together identical alignment columns
```
pos     0123456789...
seq1    ACACTCAGTC...
seq2    ACACTCAGTC...
seq3    ATACTCTGTC...
seq4    ATACTCTGTC...
```
In this example, columns 0,2 are identical, as are 3,5,9 and 4,8.
Each of these groups of identical columns would require exactly the same calculation.
Similarly, columns that are all one state up to unknown sites will never be inferred to be a state other than that one state.

If such identical alignment columns are common, substantial savings can be achieved by compressing by grouping identical columns.
If almost every column is diverse, the overhead might cost more than is gained.

### Multiple types of sequences

For flu, genomes come in segments -- should these sequences be saved as list of multiple sequences, or concatenated?
Flu genome segments all use the same alphabet, so this could be implemented either way.
But we might also want to add other discrete characters such as location data or host, which would have different alphabets.

There is thus a case for organizing sequence data into a more general structure that can hold discrete traits (sequences, locations, hosts).
In the case of sequences, these data would be quite high dimensional with an alphabet size of 4 or 5 for nucleotides, in case of hosts or locations these data would be one dimensional but with a large alphabet.

### Models
Edges in the graph are the structures along which the data change and thus along which the evolution models act.
Most of the time, we will have the same models acting across the entire tree but it is conceivable that models could change along the tree.

Each model would act either on a one-dimensional vector such as host species, or an entire matrix/profile encompassing all positions in the sequence.
In some models, some part of the sequence might evolve at different rates $\mu$ or have different equilibrium frequencies.

There is thus quite a bit of trade-off between biological realism and speed.
In most cases, fairly simple models are used.



