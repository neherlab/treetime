from collections import namedtuple
import numpy as np

# root_state is the count of each state at the root, i.e. [54, 43, 23, 12] for A, C, G, T
MutationCounts = namedtuple('MutationCounts', ['nij', 'Ti', 'root_state'])

def infer(mutation_counts_by_partition, n_states, fixed_pi=None, share_pi=False, pc=1.0, dp=1e-5, Nit=40):
    from scipy import linalg as LA

    # matrix of pseudo-counts for each state pair, added to nij to avoid zero counts
    pc_mat = pc * np.ones((n_states, n_states))
    pc_row = pc * np.ones(n_states)*(n_states-1)
    np.fill_diagonal(pc_mat, 0.0)

    W_ij_numerator = np.sum(mc.nij + mc.nij.T for mc in mutation_counts_by_partition) + 2 * pc_mat
    if share_pi: # one single pi for all partitions as a one-element list of vectors
        pi_numerator = [np.sum(mc.root_state + np.sum(mc.nij, axis=1) for mc in mutation_counts_by_partition) + pc_row]
    else: # one pi per partition as a list of vectors
        pi_numerator = [mc.root_state + np.sum(mc.nij, axis=1) + pc_row for mc in mutation_counts_by_partition]

    mu_numerator = [mc.nij.sum() + pc for mc in mutation_counts_by_partition]

    length_by_partition = [mc.root_state.sum() for mc in mutation_counts_by_partition]
    total_length = sum(length_by_partition)

    iter_count = 0
    if fixed_pi:
        pi = [np.copy(fixed_pi)]
    else:
        if share_pi:
            pi = [np.ones(n_states)]
        else:
            pi = [np.ones(n_states) for _ in mutation_counts_by_partition]

    pi = [p / p.sum() for p in pi]
    mu = [(mc.nij.sum() + pc) / (mc.Ti.sum() + pc) for mc in mutation_counts_by_partition]

    W_ij_old = np.zeros_like(W_ij_numerator)
    W_ij = np.ones_like(W_ij_numerator)

    # if pi is fixed, this will immediately converge
    while LA.norm(W_ij_old - W_ij) > dp and iter_count < Nit:
        iter_count += 1
        W_ij_old = np.copy(W_ij)

        if share_pi:
            W_ij_denominator = np.sum([muk*(np.outer(pi[0], mc.Ti) + np.outer(mc.Ti, pi[0])) for muk, mc in zip(mu, mutation_counts_by_partition)], axis=0) + 2 * pc_mat
        else:
            W_ij_denominator = np.sum([muk*(np.outer(mc.Ti, pik) + np.outer(pik, mc.Ti)) for muk, pik, mc in zip(mu, pi, mutation_counts_by_partition)], axis=0) + 2 * pc_mat
        W_ij = W_ij_numerator / W_ij_denominator
        np.fill_diagonal(W_ij, 0)

        if share_pi:
            scale_factor = avg_transition(W_ij, pi, [total_length], total_length)
        else:
            scale_factor = avg_transition(W_ij, pi, length_by_partition, total_length)

        W_ij = W_ij / scale_factor

        if share_pi:
            pi_denominator = np.sum([muk * lenk * np.dot(W_ij, mc.Ti) for lenk, muk, mc in zip(length_by_partition, mu, mutation_counts_by_partition)], axis=0)/total_length + pc_row + total_length
        else:
            pi_denominator =        [muk * np.dot(W_ij, mc.Ti) + lenk + pc_row for lenk, muk, mc in zip(length_by_partition, mu, mutation_counts_by_partition)]

        # depending on whether pi is shared or not, we either have one pi vector or one per partition
        pi = [pi_numerator[i] / pi_denominator[i] for i in range(len(pi))]
        pi = [pi[k]/pi[k].sum() for k in range(len(pi))]

        if share_pi:
            mu_denominator = [np.sum(pi[0] * (W_ij.dot(mc.Ti))) + pc for mc in mutation_counts_by_partition]
        else:
            mu_denominator = [np.sum(pik * (W_ij.dot(mc.Ti))) + pc for pik, mc in zip(pi, mutation_counts_by_partition)]

        # loop over all partitions
        mu = [mu_numerator[i] / mu_denominator[i] for i in range(len(mu))]


    gtrs = []
    for i, mc in enumerate(mutation_counts_by_partition):
        gtrs.append([mu[i], W_ij, pi[0] if share_pi else pi[i]])
    return gtrs

def avg_transition(W_ij, pi, length_by_partition, total_length):
    return np.sum([lenk * np.sum(pik * (W_ij.dot(pik))) for lenk, pik in zip(length_by_partition, pi)])/total_length




def test_infer():
    from packages.legacy.treetime.treetime import gtr
    from scipy import linalg as LA

    nij = np.array([[0, 10, 5, 2], [10, 0, 3, 1], [5, 3, 0, 4], [2, 1, 4, 0]])
    Ti = np.array([100, 100, 100, 100])
    root_state = np.array([50, 50, 50, 50])
    mutation_counts = MutationCounts(nij=nij, Ti=Ti, root_state=root_state)

    gtrs = infer([mutation_counts], n_states=4)
    print(gtrs[0].Q)

def test_multiple_partitions():
    from packages.legacy.treetime.treetime import gtr
    from scipy import linalg as LA

    nij1 = np.array([[0, 10, 5, 2], [10, 0, 3, 1], [5, 3, 0, 4], [2, 1, 4, 0]])
    Ti1 = np.array([100, 100, 100, 100])
    root_state1 = np.array([50, 50, 50, 50])
    mutation_counts1 = MutationCounts(nij=nij1, Ti=Ti1, root_state=root_state1)

    nij2 = np.array([[0, 20, 10, 4], [20, 0, 6, 2], [10, 6, 0, 8], [4, 2, 8, 0]])
    Ti2 = np.array([120, 89, 120, 87])
    root_state2 = np.array([100, 100, 100, 100])
    mutation_counts2 = MutationCounts(nij=nij2, Ti=Ti2, root_state=root_state2)

    gtrs = infer([mutation_counts1, mutation_counts2], n_states=4)
    print(gtrs[0])
    print(gtrs[1])