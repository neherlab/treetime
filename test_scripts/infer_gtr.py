import numpy as np


def avg_transition(W, pi) -> np.float64:
    return 0.0


def distance(pi_old, pi) -> np.float64:
    return sum([(x_i - y_i) ** 2 for x_i, y_i in zip(pi_old, pi)]) ** 0.5


def infer_gtr(
    nij, Ti, root_state, fixed_pi=None, alphabet=None, pc=1.0, dp=1e-5, Nit=40
):
    r"""
    Infer a GTR model by specifying the number of transitions and time spent in each
    character. The basic equation that is being solved is

    :math:`n_{ij} = pi_i W_{ij} T_j`

    where :math:`n_{ij}` are the transitions, :math:`pi_i` are the equilibrium
    state frequencies, :math:`W_{ij}` is the "substitution attempt matrix",
    while :math:`T_i` is the time on the tree spent in character state
    :math:`i`. To regularize the process, we add pseudocounts and also need
    to account for the fact that the root of the tree is in a particular
    state. the modified equation is

    :math:`n_{ij} + pc = pi_i W_{ij} (T_j+pc+root\_state)`

    Parameters
    ----------

     nij : nxn matrix
        The number of times a change in character state is observed
        between state j and i

     Ti :n vector
        The time spent in each character state

     root_state : n vector
        The number of characters in state i in the sequence
        of the root node.

     pc : float
        Pseudocounts, this determines the lower cutoff on the rate when
        no substitutions are observed

    fixed_pi : n vector of None.

    dp:  convergence criterion

    Nit: maximum number of iterations

    alphabet : Alphabet

    """
    # assert that size of n_ij is nxn and size of Ti is n, where n is the alphabet size

    pc_mat = pc * np.ones_like(nij)
    np.fill_diagonal(pc_mat, 0.0)
    np.fill_diagonal(nij, 0.0)
    pi_old = np.zeros_like(Ti)

    if fixed_pi is None:
        pi = np.ones_like(Ti)
    else:
        pi = np.copy(fixed_pi)

    pi /= pi.sum()
    W_ij = np.ones_like(nij)
    mu = (nij.sum() + pc) / (Ti.sum() + pc)  # initial guess for the rate

    # if pi is fixed, this will immediately converge
    iteration_counter = 0
    while distance(pi_old, pi) > dp and interaction_counter < Nit:
        iteration_counter += 1
        pi_old = np.copy(pi)
        W_ij = (
            (nij + nij.T + 2 * pc_mat)
            / mu
            / (np.outer(pi, Ti) + np.outer(Ti, pi) + 2 * pc_mat)
        )  # np.outer(x,y) of two vectors x and y is a matrix where the i,j element is x_i*y_j

        np.fill_diagonal(W_ij, 0)
        scale_factor = avg_transition(W_ij, pi)  # this function already exists

        W_ij = W_ij / scale_factor
        if fixed_pi is None:
            pi = (np.sum(nij + pc_mat, axis=1) + root_state) / (
                mu * np.dot(W_ij, Ti) + root_state.sum() + np.sum(pc_mat, axis=1)
            )
            pi /= pi.sum()
            mu = (nij.sum() + pc) / (np.sum(pi * (W_ij.dot(Ti))) + pc)
        else:
            mu = (nij.sum() + pc) / (np.sum(pi * (W_ij.dot(pi))) * Ti.sum() + pc)

    if count >= Nit:
        if distance(pi_old, pi) > dp:
            gtr.logger("the iterative scheme has not converged", 3, warn=True)
        elif np.abs(1 - np.max(pi.sum(axis=0))) > dp:
            gtr.logger(
                "the iterative scheme has converged, but proper normalization was not reached",
                3,
                warn=True,
            )

    gtr.assign_rates(mu=mu, W=W_ij, pi=pi)
    return gtr
