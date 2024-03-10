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

L = len(seqs['A'])

def fitch_state(child_states):
    intersection = set.intersection(*child_states)
    if intersection:
        return intersection
    else:
        return set.union(*child_states)

## postorder traversal, deal with leaves first
for n in tree.get_terminals():
    n.seq = seqs[n.name]
    n.state = [set([nuc]) for nuc in n.seq]

## postorder traversal, deal with internal nodes
for n in tree.get_nonterminals(order='postorder'):
    n.state = [fitch_state([c.state[pos] for c in n.clades]) for pos in range(L)]

## calculate the root state
from random import choice
tree.root.inferred_seq = ''.join([choice(list(n.state[pos])) for pos in range(L)])

## preorder traversal. 
for n in tree.get_nonterminals(order='preorder'):
    for c in n.clades:
        c.inferred_seq = ''.join([n.inferred_seq[pos] if n.inferred_seq[pos] in c.state[pos] else choice(list(c.state[pos])) for pos in range(L)])


# compare with treetime 

from treetime import TreeAnc
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
aln = MultipleSeqAlignment([SeqRecord(seq=seqs[k], id=k) for k in seqs])
tt = TreeAnc(tree=Phylo.read(StringIO("((A:0.1,B:0.2):0.1,(C:0.2,D:0.12):0.05):0.01;"), 'newick'), aln=aln, gtr=myGTR, compress=False)

tt.infer_ancestral_sequences(method='fitch')

for n1,n2 in zip(tt.tree.find_clades(), tree.find_clades()):
    print(tt.sequence(n1), n2.inferred_seq,tt.sequence(n1)==n2.inferred_seq)
