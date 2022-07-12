
# Sequence evolution models

Sequence evolution models describe the stochastic ways sequences change over time do to mutations.
If we consider the just a single site and denote the probability of finding that site in state $i$ as $p_i$, these models take the form
$$
\frac{d p_i}{dt} = \sum_j Q_{ij} p_j
$$
where the matrix $Q_{ij}$ describes the rate at which the mutation changes the site from state $j$ to $i$.
This is a linear equation and it is solved by
$$
p_i(t) = \left(e^{\mathbf{Q} t}\right)_{ij} p_j(0)
$$
where $\left(e^{\mathbf{Q} t}\right)_{ij}$ is the $i,j$ element of the matrix exponential of $\mathbf{Q} t$.
This matrix exponential is central to sequence evolution and one of the most computationally demanding parts.

Since the matrix $Q_{ij}$ describes the evolution of a probability distribution, it needs to conserve probability.
In other words the sum $\sum_i p_i=1$ and
$$
\sum_i \frac{d p_i}{dt} = 0 = \sum_{ij} Q_{ij} p_j
$$
Since this has to be valid for any $p_j$, we need $\sum_{i} Q_{ij} = 0$ which is typically enforced by requiring $Q_{jj} = - \sum_{i\neq j} Q_{ij}$.

This has the additional consequence that there has to be at least one vanishing eigenvalue.
Typically, this is exactly one eigenvalue but there could be special non-reducible models that have more than one (we have never considered such cases).
We adopt the convention that the vanishing eigenvalue is $\lambda_n$ where $n$ is the size of the alphabet.
The corresponding right eigenvalue $v^n_i$ is the equilibrium distribution $p_i$.
$$
\frac{d p_i}{dt} 0 = \sum_j Q_{ij} \pi_j
$$
Note that the corresponding left eigenvalue is a vector $\bar{v}^n_i$ with just $\bar{v}^n_i = 1$ at every entry.
$$
 0 = \sum_i \bar{v}^n_i Q_{ij} = \sum_i  Q_{ij} = 0
$$
The left and right eigenvectors are ortho-normal, we can further deduce that
$$
\sum_i \bar{v}^{n}_i v^m_i = \sum_i v^m_i = 0 \quad \mathrm{m\neq n}
$$
that all right eigenvectors with non-zero eigenvalues have to sum to zero.

## General time reversible models

Most models that are used in phylogenetics are so called time reversible models with the property that in equilibrium there are no net fluxes, that is that the flux from state $i$ to $j$ equals the flux from $j$ to $i$:
$$
Q_{ij} \pi_j = Q_{ji}\pi_i
$$
A simple and general way to enforce this is to parameterize the rate matrix as
$$
Q_{ij} = \pi_i W_{ij} \quad \mathrm{i\neq j}
$$
where $W_{ij}$ is a symmetric matrix.
The diagonal $Q_jj = -\sum_i Q_{ij}$ as before.

### Eigenvector calculations for general GTR models

Calculating left and right eigenvectors of a general matrix is numerically challenging.
But for GTR models, this problem can be reduced to the calculation of eigenvectors of a symmetric matrix.
Instead of calculating eigenvectors of $Q_{ij}$, we define the diagonal matrix $D_{ij} = \delta_{ij}\sqrt{\pi_j}$ and consider
$$
\tilde{\mathbf{Q}} = \mathbf{D^{-1} Q D} = \sqrt{\pi_i} W_{ij}\sqrt{\pi_j}
$$
We can now calculate eigenvalues and eigenvectors $\tilde{\mathbf{Q}}$, for which being a symmetric matrix the right and left eigenvectors agree and standard routines such as \texttt{eigh} can be used.
An eigenvector $w^k_i$ to $\tilde{\mathbf{Q}}$ obeys
$$
\mathbf{D} \mathbf{\tilde{Q}}  \mathbf{D}^{-1} \mathbf{D} w^k_i = \mathbf{D} \lambda_k  w^k_i = \mathbf{Q}  \mathbf{D} w^k_i
$$
and $\mathbf{D} w^k_i$ is thus a right eigenvector to $\mathbf{Q}$ with eigenvector $\lambda_k$.
Similarly, the left eigenvectors are $\mathbf{D}^{-1} w^k_i$.













