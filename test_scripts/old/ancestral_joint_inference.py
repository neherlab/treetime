import numpy as np
from treetime import GTR
from treetime.seq_utils import seq2prof, profile_maps, prof2seq
from Bio import Phylo
from io import StringIO

tree = Phylo.read(StringIO("((A:0.1,B:0.2):0.1,(C:0.2,D:0.12):0.05):0.01;"), 'newick')

seqs = {'A':'ACATCGCC',
        'B':'ACATCCCT',
        'C':'ACGGCCCT',
        'D':'ACGGCCCT'}

myGTR = GTR.standard('JC69', alphabet='nuc')

# Barebones implementation of marginal likelihood calculation and ancestral sequence inference
# this uses elementary operations from the GTR model but no other treetime functionality
# it treats sequences as fixed length arrays and does not compress identical columns etc
# each node keeps copies of its profile, the message from its parent, and the message to its parent
# this can use a lot of memory 


## postorder traversal, deal with leaves first
for n in tree.get_terminals():
    n.seq = seqs[n.name]
    eQt = myGTR.expQt(n.branch_length)
    n.LH_given_downstream = seq2prof(n.seq, profile_map=profile_maps['nuc'])

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
        parent_msg = np.prod([n.msg_from_parent] + [s.msg_to_parent for s in n.clades if s!=c], axis=0)
        c.msg_from_parent = myGTR.evolve(parent_msg, c.branch_length)
        c.marginal_LH = c.msg_from_parent*c.profile
        c.marginal_LH /= c.marginal_LH.sum(axis=1)[:,None]

## calculate the inferred sequences
for n in tree.find_clades():
    # convert back to sequences by picking the most likely base at each position
    # the marginal_LH are already normalized, no need to normalize again
    n.inferred_seq = ''.join(prof2seq(n.marginal_LH, myGTR, normalize=False)[0])

# compare with treetime 

from treetime import TreeAnc
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
aln = MultipleSeqAlignment([SeqRecord(seq=seqs[k], id=k) for k in seqs])
tt = TreeAnc(tree=Phylo.read(StringIO("((A:0.1,B:0.2):0.1,(C:0.2,D:0.12):0.05):0.01;"), 'newick'), aln=aln, gtr=myGTR, compress=False)

tt.infer_ancestral_sequences(marginal=True)

print("\ncompare profiles:")
print(np.abs(tt.tree.root.marginal_profile - tree.root.marginal_LH).sum())

print("\ncompare LH:")
print(tt.sequence_LH() - np.log(tree.total_LH))

for n1,n2 in zip(tt.tree.find_clades(), tree.find_clades()):
    print(tt.sequence(n1), n2.inferred_seq,tt.sequence(n1)==n2.inferred_seq)
