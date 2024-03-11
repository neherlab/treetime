import numpy as np
from treetime import GTR
from treetime.seq_utils import seq2prof, profile_maps, prof2seq
from Bio import Phylo
from io import StringIO

tree = Phylo.read(StringIO("((A:0.1,B:0.2):0.1,(C:0.2,D:0.12):0.05):0.01;"), 'newick')

seqs = {'A':'ACATCGCC',
        'B':'ACATGCCT',
        'C':'ACGGCCCT',
        'D':'GCGGCCCT'}

myGTR = GTR.standard('JC69', alphabet='nuc_nogap')

# Simple implementation of the branch length optimization without using einsum

#### SAME CODE AS IN MARGINAL INFERENCE JUST TO GENERATE THE NECESSARY PROFILES
# Only difference is that we are saving the parent_msg on the node for convenience

## postorder traversal, deal with leaves first
for n in tree.get_terminals():
    n.seq = seqs[n.name]
    n.profile = seq2prof(n.seq, profile_map=profile_maps['nuc_nogap'])
    n.msg_to_parent = myGTR.propagate_profile(n.profile, n.branch_length)

## postorder traversal, deal with internal nodes
for n in tree.get_nonterminals(order='postorder'):
    n.profile = np.prod([c.msg_to_parent for c in n.clades], axis=0)
    n.msg_to_parent = myGTR.propagate_profile(n.profile, n.branch_length)

## calculate the root state
tree.root.msg_from_parent = np.array([myGTR.Pi for i in range(len(tree.root.profile))])
tree.root.marginal_LH = tree.root.msg_from_parent*tree.root.profile

# if not careful, this will result in underflows for long sequences
tree.total_LH = tree.root.marginal_LH.sum(axis=1).prod()

# this is a profile normalization step that we need to do a lot
tree.root.marginal_LH /= tree.root.marginal_LH.sum(axis=1)[:,None]

## preorder traversal. 
for n in tree.get_nonterminals(order='preorder'):
    for c in n.clades:
        c.parent_msg = np.prod([n.msg_from_parent] + [s.msg_to_parent for s in n.clades if s!=c], axis=0)
        c.msg_from_parent = myGTR.evolve(c.parent_msg, c.branch_length)
        c.marginal_LH = c.msg_from_parent*c.profile
        c.marginal_LH /= c.marginal_LH.sum(axis=1)[:,None]

#############################

## set up treetime to compare
from treetime import TreeAnc
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
aln = MultipleSeqAlignment([SeqRecord(seq=seqs[k], id=k) for k in seqs])
tt = TreeAnc(tree=Phylo.read(StringIO("((A:0.1,B:0.2):0.1,(C:0.2,D:0.12):0.05):0.01;"), 'newick'), aln=aln, gtr=myGTR, compress=False, alphabet='nuc_nogap')
tt.infer_ancestral_sequences(marginal=True)
#############


child_node = tree.get_terminals()[1]
# multiplicity is needed by python treetime to compress identical columns in the alignment, we will deal with this differently later in rust. 
multiplicity=np.ones(len(child_node.profile))

# implementation of prob_t 
def prob_t(parent, child, t, return_log=True, ignore_gaps=False):
    eQt = myGTR.expQt(t) # this function we already have
    res = [c.dot(eQt).dot(p) for p,c in zip(parent, child)]  # this is what we previously did with einsum
    # TODO: ignore gaps -- probably good to keep as argument to function, but need to think through implementation
    if ignore_gaps:
        print('not implemented')
    if return_log:
        return np.sum(np.log(res))
    else:
        return np.prod(res)

# test for prob_t
print("length, treetime, prob_t")
for t in [0, 0.1, 0.2, 0.3]:
    print(t, myGTR.prob_t_profiles((child_node.parent_msg, child_node.profile), multiplicity, t, return_log=True, ignore_gaps=False), 
          prob_t(child_node.parent_msg, child_node.profile, t))

# finding the t that maximixes prob_t

MAX_BRANCH_LENGTH = 10  # magic number

# cost function of square root of branch length to avoid having to deal with bounds
def cost(s, parent, child):
    return -prob_t(parent, child, s**2, return_log=True, ignore_gaps=False)

# l0 is the starting length
def optimal_t(parent, child, l0, multiplicity):
    from scipy.optimize import minimize_scalar

    res = minimize_scalar(cost, bracket=[-np.sqrt(MAX_BRANCH_LENGTH), np.sqrt(l0), np.sqrt(MAX_BRANCH_LENGTH)],
                          args=(parent, child), method='brent', tol=1e-6)
    if res['success']:
        return res['x']**2
    else:
        return np.nan
print('\n')
print("optimal_t:", optimal_t(child_node.parent_msg, child_node.profile, 0.1, multiplicity))
print("TreeTime", myGTR.optimal_t_compressed((child_node.parent_msg, child_node.profile), multiplicity=np.ones(len(child_node.profile)), profiles=True))

