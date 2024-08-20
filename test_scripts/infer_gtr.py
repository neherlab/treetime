import numpy as np
from Bio import AlignIO
import numpy as np
from lib import Graph, graph_from_nwk_str, GraphNodeBackward, GraphNodeForward
from lib import SparseSeqDis, VarPos
from treetime import GTR
from payload import NodePayload, EdgePayload
from profile_map import profile_map
from fitch import init_sequences_sparse


def avg_transition(W, pi) -> np.float64:
    return pi.dot(W.dot(pi))


def distance(pi_old, pi) -> np.float64:
    #return sum([(x_i - y_i) ** 2 for x_i, y_i in zip(pi_old, pi)]) ** 0.5
    return np.sum((pi_old-pi)**2)**0.5  # numpy version ;)


def get_mutation_counts(graph):
    root = graph.get_one_root().payload()
    nij = np.zeros((4, 4))
    Ti = np.zeros(4)
    seq = root.sparse_sequences[0].seq
    root_state = np.array([seq.composition[nuc] for nuc in 'ACGT'])
    char_to_index = {nuc: i for i, nuc in enumerate('ACGT')}
    for e in graph.get_edges():
        target_seq = graph.get_node(e.target()).payload().sparse_sequences[0].seq
        bl = e.payload().branch_length
        for i, nuc in enumerate('ACGT'):
            Ti[i] += bl * target_seq.composition[nuc]

        for mut in e.payload().sparse_sequences[0].muts:
            i, j = char_to_index[mut.qry], char_to_index[mut.ref]
            nij[i,j] += 1
            Ti[i] -= 0.5*bl
            Ti[j] += 0.5*bl

    return nij, Ti, root_state


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
    while distance(pi_old, pi) > dp and iteration_counter < Nit:
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

    if iteration_counter >= Nit:
        if distance(pi_old, pi) > dp:
            print("the iterative scheme has not converged")
        elif np.abs(1 - np.max(pi.sum(axis=0))) > dp:
            print("the iterative scheme has converged, but proper normalization was not reached")

    return {"W": W_ij, "pi": pi, "mu": mu}


def test():
    nij = np.array([[0, 1, 2, 1], [1, 0, 3, 2], [2, 3, 0, 1], [2, 3, 3, 0]])
    Ti = np.array([12.0, 20., 14.0, 12.4])
    root_state = np.array([3, 2, 3, 4])
    print(infer_gtr(nij, Ti, root_state, pc=0.1))

if __name__=="__main__":
    fname_nwk = 'data/ebola/ebola.nwk'
    fname_seq = 'data/ebola/ebola_dna.fasta'
    fname_nwk = 'test_scripts/data/tree.nwk'
    fname_seq = 'test_scripts/data/sequences.fasta'
    with open(fname_nwk) as fh:
        nwkstr = fh.read()
    G = graph_from_nwk_str(nwk_string=nwkstr, node_payload_factory=NodePayload, edge_payload_factory=EdgePayload)

    aln = {seq.id: str(seq.seq).upper() for seq in AlignIO.read(fname_seq, 'fasta')}
    gtr = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')
    init_sequences_sparse(G, [aln], [gtr])


    nij, Ti, root_state = get_mutation_counts(G)

    print(infer_gtr(nij, Ti, root_state, pc=0.1))
