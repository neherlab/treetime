import numpy as np 


# sequence representation: for large trees, we should avoid storing the entire sequence at every node. 
# since along every branch, only a few positions change, there should be efficient ways to handle this

# in addition, many sequences are incomplete, so we should record ranges of missing information


# from an evolutionary perspective, different parts of the sequence might evolve via different models
# this is typically achieved by partitioning the sequence into different regions each with their separate evolutionary models

# here, we sketch the compressed representation for a single partition. 


from Bio import Phylo ,SeqIO
from io import StringIO

tree_fname = '../test/treetime_examples/data/ebola/ebola.nwk'
aln_fname = '../test/treetime_examples/data/ebola/ebola.fasta'

dummy=True

if dummy:
    tree = Phylo.read(StringIO("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;"), 'newick')
    # in practice, we don't want all sequences in memory. Idealy we'd stream then in the same order as leaves of the tree. 
    seqs = {'A':'ACATCGCCNNA--G',
            'B':'GCATCCCTGTA-TG',
            'C':'CCGGCGATGTATTG',
            'D':'TCGGCCGTGTRTTG'}
else:
    # Ebola test data
    tree = Phylo.read(tree_fname, 'newick')
    # remove U in favor or T
    seqs = {r.id:str(r.seq.upper().replace('U', 'T')) for r in SeqIO.parse(aln_fname, 'fasta')}

# pull out the sequence length for convenience
L = len(seqs[tree.get_terminals()[0].name])

# utility functions to determine the ranges of missing information and gaps
def find_ambiguous_ranges(seq):
    return find_char_ranges(seq, 'N')

def find_gap_ranges(seq):
    return find_char_ranges(seq, '-')

def find_char_ranges(seq, char):
    ranges = []
    start = None
    for pos, nuc in enumerate(seq):
        if nuc != char and (start is not None):
            ranges.append((start, pos))
            start = None
        elif nuc == char and start is None:
            start = pos
    return ranges

def range_intersection(range_sets):
    if any([len(r)==0 for r in range_sets]):
        return []
    current_ranges = list(range_sets[0])
    for next_ranges in range_sets[1:]:
        new_ranges = []
        ri1, ri2 = 0, 0
        r1 = current_ranges[ri1]
        r2 = next_ranges[ri2]
        while ri1<len(current_ranges) and ri2<len(next_ranges):
            if r2[0]>r1[1]:
                ri1 += 1
                if ri1<len(current_ranges):
                    r1 = current_ranges[ri1]
            elif r1[0]>r2[1]:
                ri2 += 1
                if ri2<len(next_ranges):
                    r2 = next_ranges[ri2]
            else:
                new_ranges.append((max(r1[0], r2[0]), min(r1[1], r2[1])))
                if r1[1]<r2[1]:
                    ri1 += 1
                    if ri1<len(current_ranges):
                        r1 = current_ranges[ri1]
                else:
                    ri2 += 1
                    if ri2<len(next_ranges):
                        r2 = next_ranges[ri2]
        current_ranges = new_ranges
    return current_ranges

def ranges_contain(ranges,pos):
    return any([r[0]<=pos<r[1] for r in ranges])

# definition and function to handle mixed sites
mixed_sites = {'R':set('AG'), 'Y':set('CT'), 
               'S':set('GC'), 'W':set('AT'), 'K':set('GT'), 'M':set('AC'), 
               'B':set('CGT'), 'D':set('AGT'), 'H':set('ACT'), 'V':set('ACG')}

def find_mixed_sites(seq):
    return {pos:nuc for pos, nuc in enumerate(seq) if nuc not in 'ACGTN-'}, {pos:mixed_sites[nuc] for pos, nuc in enumerate(seq) if nuc not in 'ACGTN-'}

# Algorithm to compactly represent the sequence at each node via mutations on the tree while avoiding 
# storing all sequences at the same time (this would only work if we can stream them in the right order, 
# but this is also a warm up exercise for the an approximate marginal ancestral inference). 

## postorder traversal
for n in tree.find_clades(order='postorder'):
    if n.is_terminal():
        # at each terminal node, temporarily store the sequence and ranges of N, - and mixed sites
        n.seq = seqs[n.name]
        n.ambiguous = find_ambiguous_ranges(n.seq)
        n.gaps = find_gap_ranges(n.seq)
        # n.mixed stores the exact character at each mixed positions, the non_consensus stores the possible states
        n.mixed, n.non_consensus = find_mixed_sites(n.seq)
    else:
        # positions that are N or - in all children are still N or - in the parent
        n.ambiguous = range_intersection([c.ambiguous for c in n.clades])
        n.gaps = range_intersection([c.gaps for c in n.clades])
        # all sites that are not N or - but not fixed will need special treatment
        non_consensus_positions = set.union(*[set(c.non_consensus.keys()) for c in n.clades])
        n.non_consensus = {}
        n.seq = ['']*L # construct sequence of node, will be deleted later again
        for pos in range(L):
            # skip ambiguous and gaps
            if ranges_contain(n.ambiguous, pos):
                continue
            if ranges_contain(n.gaps, pos):
                continue
            # deal with positions that are variable in at least one child
            if pos in non_consensus_positions:
                # indeterminate in at least one child
                isect = set.intersection(*[c.non_consensus.get(pos, set({c.seq[pos]})) 
                                           for c in n.clades if not (ranges_contain(c.ambiguous, pos) or ranges_contain(c.gaps, pos))])
                
                if len(isect)==1:
                    n.seq[pos] = isect.pop()
                elif len(isect)>1:
                    n.non_consensus[pos] = isect
                    for c in n.clades:
                        if pos not in c.non_consensus:
                            c.non_consensus[pos] = set([c.seq[pos]])
                else:
                    n.non_consensus[pos] = set.union(*[c.non_consensus.get(pos, set([c.seq[pos]])) 
                                           for c in n.clades if not (ranges_contain(c.ambiguous, pos) or ranges_contain(c.gaps, pos))])
                    for c in n.clades:
                        if pos not in c.non_consensus:
                            c.non_consensus[pos] = set([c.seq[pos]])
            else: # deal with all other positions. 
                # this could probably be sped up by explicitly checking whether the states of all children are equal. 
                states = set([c.seq[pos] for c in n.clades])
                if len(states)==1: # if all children are equal
                    n.seq[pos] = states.pop()
                else: # if children differ
                    n.non_consensus[pos] = states
                    for c in n.clades:
                        c.non_consensus[pos] = set([c.seq[pos]])
        # no longer need sequences of children
        for c in n.clades:
            del c.seq

# determine the sequence at the root
for pos, states in tree.root.non_consensus.items():
    tree.root.seq[pos] = states.pop()  # should be random choice


# we now have a complete sequence at the root and should be able to delete the non_consensus 
# there might still be ambiguous positions, but that means we have no information anywhere...

tree.root.tmp_seq = list(tree.root.seq)
tree.root.mutations = {}

# do a pre-order traversal to construct all necessary mutations
for n in tree.get_nonterminals(order='preorder'):
    for c in n.clades:
        c.mutations = {}
        # we need this temporary sequence only for internal nodes
        if not c.is_terminal():
            c.tmp_seq = list(n.tmp_seq)

        # all positions that potentially differ in the child c from the parent n are in `c.non_consensus` 
        for pos, states in c.non_consensus.items():
            if n.tmp_seq[pos] not in states:
                # in this case we need a mutation to one state in states
                state = states.pop()
                c.mutations[pos] = (n.tmp_seq[pos], state)
                if not c.is_terminal():
                    c.tmp_seq[pos] = state
    # no longer needed
    del n.non_consensus
    del n.tmp_seq

# function to reconstruct the sequence of a node from the sequence at the root and mutations
def reconstruct_seq(tree, node):
    # get the path in the tree from the root to node
    path = tree.get_path(n)
    # copy the root sequence as array to be able to modify by element and slice
    seq = np.array(tree.root.seq)
    # apply mutations
    for anc in path:
        for pos in anc.mutations:
            seq[pos] = anc.mutations[pos][1]
    # introduce N, gaps, and mixed sites
    for r in node.ambiguous:
        seq[r[0]:r[1]] = 'N'
    for r in node.gaps:
        seq[r[0]:r[1]] = '-'
    for pos, state in n.mixed.items():
        seq[pos] = state
    return ''.join(seq)

fails = []
for n in tree.get_terminals():
    seq = reconstruct_seq(tree, n)
    if seq != seqs[n.name]:
        fails.append(n.name)
        print(n.name, seqs[n.name]==seq)

if len(fails):
    print("DAMN: the following sequence were not reconstructed correctly:", fails)
else:
    print("TADA: all sequences were reconstructed correctly.")




