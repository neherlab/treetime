import numpy as np 
from collections import defaultdict

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
    nwk_str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;"
    tree = Phylo.read(StringIO(nwk_str), 'newick')
    # in practice, we don't want all sequences in memory. Idealy we'd stream then in the same order as leaves of the tree. 
    seqs = {'A':'ACATCGCCNNA--G',
            'B':'GCATCCCTGTA-NG',
            'C':'CCGGCGATGTATTG',
            'D':'TCGGCCGTGTRTTG'}
    
    # seqs = {'A':'ACATCGCCTTATTGAGGT',
    #         'B':'GCATCCCTGTATTGAGGT',
    #         'C':'CCGGCGATGTATTGAGGT',
    #         'D':'TCGGCCGTGTTTTGAGGT'}
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

def find_undetermined_ranges(seq):
    # this should be the union of the gap and N ranges
    return find_char_ranges(seq.replace('-', 'N'), 'N')

def find_char_ranges(seq, char):
    ranges = []
    start = None
    for pos, nuc in enumerate(seq):
        if nuc != char and (start is not None):
            ranges.append((start, pos))
            start = None
        elif nuc == char and start is None:
            start = pos
    if start:
        ranges.append((start, len(seq)))
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
        n.undetermined = find_undetermined_ranges(n.seq)
        n.gaps = find_gap_ranges(n.seq)
        # n.mixed stores the exact character at each mixed positions, the non_consensus stores the possible states
        n.mixed, n.non_consensus = find_mixed_sites(n.seq)
    else:
        # positions that are N or - in all children are still N or - in the parent
        n.undetermined = range_intersection([c.undetermined for c in n.clades])
        # all sites that are not N or - but not fixed will need special treatment
        non_consensus_positions = set.union(*[set(c.non_consensus.keys()) for c in n.clades])
        n.non_consensus = {}
        n.seq = ['']*L # construct sequence of node, will be deleted later again
        for pos in range(L):
            # skip ambiguous and gaps
            if ranges_contain(n.undetermined, pos):
                continue
            # deal with positions that are variable in at least one child
            if pos in non_consensus_positions:
                # indeterminate in at least one child
                isect = set.intersection(*[c.non_consensus.get(pos, set({c.seq[pos]})) 
                                           for c in n.clades if not ranges_contain(c.undetermined, pos)])
                if '-' in isect:
                    import ipdb; ipdb.set_trace()
                if len(isect)==1:
                    n.seq[pos] = isect.pop()
                elif len(isect)>1:
                    n.non_consensus[pos] = isect
                    for c in n.clades:
                        if pos not in c.non_consensus:
                            c.non_consensus[pos] = set([c.seq[pos]])
                else:
                    n.non_consensus[pos] = set.union(*[c.non_consensus.get(pos, set([c.seq[pos]])) 
                                           for c in n.clades if not ranges_contain(c.undetermined, pos)])
                    for c in n.clades:
                        if pos not in c.non_consensus:
                            c.non_consensus[pos] = set([c.seq[pos]])
            else: # deal with all other positions. 
                # this could probably be sped up by explicitly checking whether the states of all children are equal. 
                states = set([c.seq[pos] for c in n.clades if not ranges_contain(c.undetermined, pos)])
                if len(states)==1: # if all children are equal
                    n.seq[pos] = states.pop()
                else: # if children differ
                    n.non_consensus[pos] = states
                    for c in n.clades:
                        c.non_consensus[pos] = set([c.seq[pos]])
        print(len(n.non_consensus))
        # no longer need sequences of children
        for c in n.clades:
            del c.seq

# def dump_tree(tree):
#     node_dicts = []
#     for n in tree.find_clades(order='postorder'):
#         node_dicts.append({myGTR.expQt
#             "name": n.name,
#             "is_terminal": n.is_terminal(),
#             "mutations": n.__dict__.get("mutations"),
#             "gaps": n.__dict__.get("gaps"),
#             "ambiguous": n.__dict__.get("ambiguous"),
#             "mixed": n.__dict__.get("mixed"),
#             "non_consensus": list(map(lambda x: (x[0], list(x[1])), n.__dict__.get("non_consensus").items())),
#         })
#     import json
#     return json.dumps(node_dicts, indent=2)
#
# print(dump_tree(tree))

# determine the sequence at the root
for pos, states in tree.root.non_consensus.items():
    tree.root.seq[pos] = states.pop()  # should be random choice


# we now have a complete sequence at the root and should be able to delete the non_consensus 
# there might still be ambiguous positions, but that means we have no information anywhere...

tree.root.tmp_seq = list(tree.root.seq)
tree.root.muts = {}

# do a pre-order traversal to construct all necessary mutations
for n in tree.get_nonterminals(order='preorder'):
    for c in n.clades:
        c.muts = {}
        # we need this temporary sequence only for internal nodes
        if not c.is_terminal():
            c.tmp_seq = list(n.tmp_seq)

        # all positions that potentially differ in the child c from the parent n are in `c.non_consensus` 
        for pos, states in c.non_consensus.items():
            if n.tmp_seq[pos] not in states:
                # in this case we need a mutation to one state in states
                state = states.pop()
                c.muts[pos] = (n.tmp_seq[pos], state)
                if not c.is_terminal():
                    c.tmp_seq[pos] = state
    # no longer needed
    del n.non_consensus
    del n.tmp_seq


# function to reconstruct the sequence of a node from the sequence at the root and mutations
def reconstruct_raw_seq(tree, node):
    # get the path in the tree from the root to node
    path = tree.get_path(node)
    # copy the root sequence as array to be able to modify by element and slice
    seq = np.array(tree.root.seq)
    # apply mutations
    for anc in path:
        for pos in anc.muts:
            seq[pos] = anc.muts[pos][1]
    return seq

# function to reconstruct the sequence of a node from the sequence at the root and mutations
def reconstruct_seq(tree, node):
    seq = reconstruct_raw_seq(tree, node)
    introduce_non_nucs(node, seq)
    return ''.join(seq)

def introduce_non_nucs(node, seq):
    # introduce N, gaps, and mixed sites
    for r in node.ambiguous:
        seq[r[0]:r[1]] = 'N'
    for r in node.gaps:
        seq[r[0]:r[1]] = '-'
    for pos, state in node.mixed.items():
        seq[pos] = state


# iterate over sequence while keeping track of the sequence
# straight forward in pre-order or post-order traversal
# unsure how one would best do this in breadth-first. 
seq = np.array(tree.root.seq) # to allow element wise assignment

def pre_order(tree, node, seq, do_stuff):
    # pick up sequence from parent, apply mutations
    for pos, (anc, der) in node.muts.items():
        seq[pos] = der
    # execute what is necessary (this would be the yield in an iterator)
    do_stuff(tree, node, seq)
    # call the same function for the children
    for c in node.clades:
        pre_order(tree, c, seq, do_stuff)

    # undo the mutations
    for pos, (anc, der) in node.muts.items():
        seq[pos] = anc

def post_order(tree, node, seq, do_stuff):
    # pick up sequence from parent, apply mutations
    for pos, (anc, der) in node.muts.items():
        seq[pos] = der
    # execute what is necessary (this would be the yield in an iterator)
    # call the same function for the children
    for c in node.clades:
        post_order(tree, c, seq, do_stuff)
    do_stuff(tree, node, seq)

    # undo the mutations
    for pos, (anc, der) in node.muts.items():
        seq[pos] = anc

def check_seq(tree, node, seq):
    if node.is_terminal():
        seq_copy = np.copy(seq)
        introduce_non_nucs(node, seq_copy)
        seq_str=''.join(seq_copy)
        print(f"{node.name}\t", "".join(seq), seq_str, seqs[node.name]==seq_str)
    else:
        print(f"{node.name}\t", "".join(seq))

if dummy:
    pre_order(tree, tree.root, seq, check_seq)
    post_order(tree, tree.root, seq, check_seq)

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

## calculate nucleotide composition -- only for dev purposes
def nuc_comp(node, seq):
    node.nuc_composition = defaultdict(int)
    if node.is_terminal():
        seq_copy = np.copy(seq)
        introduce_non_nucs(node, seq_copy)
    else:
        seq_copy = seq

    for n in 'ACGTN-':
        node.nuc_composition[n] = np.sum(n==seq)

for n in tree.find_clades():
    seq = reconstruct_raw_seq(tree, n)
    nuc_comp(n, seq)


# with the sequence representation in place, we can now calculate the likelihood. 
from treetime import GTR
#myGTR = GTR.standard('JC69', alphabet='nuc_nogap')
myGTR = GTR.custom(pi=[0.2, 0.3, 0.15, 0.45], alphabet='nuc_nogap')
from treetime.seq_utils import profile_maps
prof_nuc = profile_maps['nuc_nogap']

# threshold to keep position among the variable ones
eps = 1e-6

# payload function that calculates the likelihood
def calc_likelihood(tree, node, seq):
    # GTR matrix associated with this branch length. Using 0 length for the root saves an extra calculation below
    node.expQt = myGTR.expQt(0 if node==tree.root else node.branch_length).T # might make sense to save this on the edge

    # we have calculated the total nucleotide composition in the sequence representation. 
    # from this, we will subtract positions that are tracked as variable positions -- hence need to copu
    fixed_nuc_count = {k:v for k,v in node.nuc_composition.items()}

    # each node will get a vector with the probability distribution of non-variable positions
    # this vector should always be peaked around the focal nucleotide and quantifies the uncertainty around it
    node.subtree_profile_fixed = {}

    if node.is_terminal():
        # For terminal nodes, we deem mixed sites and sites that have mutations in the branch to the parent as variable
        node.subtree_profile_variable = {}
        for pos, state in node.mixed.items():
            node.subtree_profile_variable[pos] = prof_nuc[state]
            fixed_nuc_count[seq[pos]] -= 1
        for pos, (anc, der) in node.muts.items():
            node.subtree_profile_variable[pos] = prof_nuc[der]
            fixed_nuc_count[seq[pos]] -= 1

        # this could be done more efficiently. We just need to look-up these positions, no need to save the flat vector.
        for rg in node.undetermined:
            for pos in range(*rg):
                node.subtree_profile_variable[pos] = prof_nuc['N']

        # node.message_to_parent = {pos: expQt.dot(prof) for pos, prof in node.subtree_profile_variable.items()}
        for ni, n in enumerate('ACGT'):
            node.subtree_profile_fixed[n] = prof_nuc[n]
            # node.subtree_profile_fixed[n] = expQt[ni, :]
    else:
        # For internal nodes, we consider all positions that are variable in any of the children
        node.subtree_profile_variable = {}
        variable_pos = set.union(*[set(c.subtree_profile_variable.keys()) for c in node.clades])
        for pos in variable_pos:
            nuc = seq[pos]
            # calculate the product of child messages
            tmp_msg = []
            for c in node.clades:
                tmp_msg.append(c.expQt.dot(c.subtree_profile_variable.get(pos, c.subtree_profile_fixed[nuc])))

            # could do selection here right away, keep only variable ones
            vec = np.prod(tmp_msg, axis=0)
            vec_norm = vec.sum()
            tree.logLH += np.log(vec_norm)

            # add position to variable states if the subleading states have a probability exceeding eps
            if vec.max()<(1-eps)*vec_norm or nuc != 'ACGT'[vec.argmax()]:
                node.subtree_profile_variable[pos] = vec/vec_norm

            # this position is accounted for, hence we can subtract it from the count of fixed nucs 
            # unless nuc is `N` or `-` since these are in nuc-composition
            if nuc and nuc in 'ACGT': 
                fixed_nuc_count[nuc] -= 1

        # collect contribution from the inert sites
        for n in 'ACGT':
            vec = np.prod([c.expQt.dot(c.subtree_profile_fixed[n]) for c in node], axis=0)
            vec_norm = vec.sum()
            tree.logLH += fixed_nuc_count[n]*np.log(vec_norm)
            node.subtree_profile_fixed[n] = vec/vec_norm

        # add position that mutate towards the parent
        for pos, (anc, der) in node.muts.items():
            if pos in node.subtree_profile_variable: 
                continue
            node.subtree_profile_variable[pos] = node.subtree_profile_fixed[der]
    # NOTE: we could save c.expQT.dot(xxx) on the edges. that would save some computation. 


# run the likelihood calculation
tree.logLH=0
post_order(tree, tree.root, tree.root.seq, calc_likelihood)

# multiply the `message_to_parent`` at the root with the equilibrium probabilties
tree.root.profile_variable = {}
inert_nucs = {k:v for k,v in tree.root.nuc_composition.items()}
# variable positions
for pos, vec in tree.root.subtree_profile_variable.items():
    tree.root.profile_variable[pos] = vec*myGTR.Pi
    vec_norm = np.sum(tree.root.profile_variable[pos])
    tree.root.profile_variable[pos]/=vec_norm
    tree.logLH += np.log(vec_norm)
    nuc = tree.root.seq[pos]
    if nuc and nuc in 'ACGT':
        inert_nucs[nuc] -= 1

# inert positions
tree.root.profile_fixed = {}
for n in 'ACGT':
    tree.root.profile_fixed[n] = tree.root.subtree_profile_fixed[n]*myGTR.Pi
    vec_norm = tree.root.profile_fixed[n].sum()
    tree.root.profile_fixed[n]/=vec_norm
    tree.logLH += inert_nucs[n]*np.log(vec_norm)


# compare with treetime
from treetime import TreeAnc
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
aln = MultipleSeqAlignment([SeqRecord(seq=Seq(seqs[k]), id=k) for k in seqs])

if dummy:
    new_tree = Phylo.read(StringIO(nwk_str), 'newick')
else:
    # Ebola test data
    new_tree = Phylo.read(tree_fname, 'newick')

tt = TreeAnc(tree=new_tree, aln=aln, gtr=myGTR, compress=False)

tt.infer_ancestral_sequences(marginal=True)

# check numerical values of the profiles at the root
eps2=0.000001
for pos in tree.root.profile_variable:
    agree = np.abs(tree.root.profile_variable[pos] - tt.tree.root.marginal_profile[pos]).sum()<eps2
    if not agree:
        print(pos, tree.root.profile_variable[pos], tt.tree.root.marginal_profile[pos])

for pos in range(L):
    if pos in tree.root.profile_variable: continue
    agree = np.abs(tree.root.profile_fixed[tree.root.seq[pos]] - tt.tree.root.marginal_profile[pos]).sum()<eps2
    if not agree:
        print(pos, tree.root.seq[pos], tree.root.profile_fixed[tree.root.seq[pos]], tt.tree.root.marginal_profile[pos])

# compare total likelihood
print(tree.logLH, tt.sequence_LH(), tree.logLH-tt.sequence_LH(), )