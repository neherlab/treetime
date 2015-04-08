
import numpy as np
from Bio import AlignIO
from time_tree import tree_anc

if __name__ == '__main__':

    alphabet = np.array(['A', 'C', 'G', 'T'])


    t = tree_anc.TreeAnc.from_file('../data/PR.B.100.nwk', 'newick')
    aln = AlignIO.read('../data/PR.B.100.fasta', 'fasta')
    t.set_seqs_to_leaves(aln)
    #t._fitch_anc()

    gtr = tree_anc.GTR(alphabet)

    # construct Jukes-Cantor GTR matrix
    # FIXME transfer to GTR constructor
    gtr.W = np.ones((alphabet.shape[0], alphabet.shape[0]))
    np.fill_diagonal(gtr.W, - ((gtr.W).sum(0) - 1))
    gtr.Pi = np.zeros(gtr.W.shape)
    np.fill_diagonal(gtr.Pi, 0.25)
    sqrtPi = np.sqrt(gtr.Pi)
    sqrtPi_inv = np.linalg.inv(sqrtPi)
    W = (sqrtPi.dot(((gtr.Pi).dot(gtr.W)))).dot(sqrtPi_inv)
    eigvals, eigvecs = np.linalg.eig(W)
    gtr.v = sqrtPi.dot(eigvecs)
    gtr.v_inv = np.linalg.inv(gtr.v)
    gtr.eigenmat = np.diagflat(eigvals)

