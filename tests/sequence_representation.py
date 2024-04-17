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
    
    # seqs = {'A':'ACATCGCCTTATAGAGGT',
    #         'B':'GCATCCCTGTATTGAGAT',
    #         'C':'CCGGCGATGTATTGAGGT',
    #         'D':'TCGGCCGTGTTTTGAGGT'}
else:
    # Ebola test data
    tree = Phylo.read(tree_fname, 'newick')
    # remove U in favor or T
    seqs = {r.id:str(r.seq.upper().replace('U', 'T')) for r in SeqIO.parse(aln_fname, 'fasta')}

seq_store = {k:np.array(list(v)) for k,v in seqs.items()}

# pull out the sequence length for convenience
L = len(seqs[tree.get_terminals()[0].name])
for n in tree.get_nonterminals():
    for c in n:
        c.parent = n
tree.root.parent=None

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
        seq = "".join(seq_store[n.name])
        n.ambiguous = find_ambiguous_ranges(seq)
        n.undetermined = find_undetermined_ranges(seq)
        n.gaps = find_gap_ranges(seq)
        # n.mixed stores the exact character at each mixed positions, the non_consensus stores the possible states
        n.mixed, n.non_consensus = find_mixed_sites(seq)
    else:
        # positions that are N or - in all children are still N or - in the parent
        n.undetermined = range_intersection([c.undetermined for c in n.clades])
        # all sites that are not N or - but not fixed will need special treatment
        non_consensus_positions = set.union(*[set(c.non_consensus.keys()) for c in n.clades])
        n.non_consensus = {}
        seq = ['?']*L # construct sequence of node, will be deleted later again

        # introduce Ns to mark indeterminate positions
        for rg in n.undetermined:
            for pos in range(*rg):
                seq[pos] = 'N'

        # indeterminate in at least one child
        for pos in non_consensus_positions:
            child_sets = []
            for c in n.clades:
                cseq = seq_store[c.name]
                if pos in c.non_consensus:
                    child_sets.append(c.non_consensus[pos])
                elif cseq[pos] in 'ACGT': # memorize child state to assign mutations
                    child_sets.append(set({cseq[pos]}))
            isect = set.intersection(*child_sets) if child_sets else set()

            if len(isect)==1:
                seq[pos] = isect.pop()
            else:
                if len(isect)>1:
                    n.non_consensus[pos] = isect
                else:
                    n.non_consensus[pos] = set.union(*child_sets)
                seq[pos] = '~'

        for pos, (nuc, child_states) in enumerate(zip(seq, zip(*[seq_store[c.name] for c in n.clades]))):
            if nuc!='?': # these positions have been dealt with above
                continue

            # this could probably be sped up by explicitly checking whether the states of all children are equal. 
            states = set([x for x in child_states if x in 'ACGT'])
            if len(states)==1: # if all children are equal
                seq[pos] = states.pop()
            else: # if children differ
                n.non_consensus[pos] = states
                seq[pos] = '~'
        seq_store[n.name] = seq

# determine the sequence at the root
seq = seq_store['root']
for pos, states in tree.root.non_consensus.items():
    seq[pos] = states.pop()  # should be random choice
# we now have a complete sequence at the root and should be able to delete the non_consensus 
# there might still be ambiguous positions, but that means we have no information anywhere...

tree.root.muts = {}
tree.root.seq = seq #this one is to keep

# do a pre-order traversal to construct all necessary mutations
for n in tree.find_clades(order='preorder'):
    if n==tree.root: continue

    n.muts = {}
    parent_seq = seq_store[n.parent.name] # read only
    seq = seq_store[n.name] # read/write
    # all positions that potentially differ in the child c from the parent n are in `c.non_consensus` 
    for pos in set.union(set(n.non_consensus.keys()), set(n.parent.non_consensus.keys())):
        if pos in n.non_consensus and parent_seq[pos] not in states:
            # in this case we need a mutation to one state in states
            state = n.non_consensus[pos].pop()
            n.muts[pos] = (parent_seq[pos], state)
            if not n.is_terminal():
                seq[pos] = state
        elif pos in n.non_consensus and parent_seq[pos] in states:
            # in this case, we can just copy the state
            seq[pos] = parent_seq[pos]
        elif parent_seq[pos]!=seq[pos]:
            n.muts[pos] = (parent_seq[pos], seq[pos])


# function to reconstruct the sequence of a node from the sequence at the root and mutations
def reconstruct_raw_seq(tree, node):
    # get the path in the tree from the root to node
    path = tree.get_path(node)
    # copy the root sequence as array to be able to modify by element and slice
    seq = np.array(tree.root.seq) # to allow element wise assignment
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
myGTR = GTR.custom(pi=[0.2, 0.3, 0.15, 0.35], alphabet='nuc_nogap')
from treetime.seq_utils import profile_maps
prof_nuc = profile_maps['nuc_nogap']

# threshold to keep position among the variable ones
eps = 1e-6

# payload function that calculates the likelihood
def subtree_profiles(tree, node, seq):
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

def calculate_root_state(tree):
    # multiply the `message_to_parent` at the root with the equilibrium probabilties
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

    # fixed positions
    tree.root.profile_fixed = {}
    for n in 'ACGT':
        tree.root.profile_fixed[n] = tree.root.subtree_profile_fixed[n]*myGTR.Pi
        vec_norm = tree.root.profile_fixed[n].sum()
        tree.root.profile_fixed[n]/=vec_norm
        tree.logLH += inert_nucs[n]*np.log(vec_norm)


def outgroup_profiles(tree, node, seq):
    if node == tree.root: # root state gets different treatment
        calculate_root_state(tree)
    else: # for other nodes, calculate outgroup_profile and profile for non-root nodes
        variable_pos = set.union(set(node.subtree_profile_variable.keys()), set(node.parent.profile_variable.keys()))
        node.outgroup_profile_variable = {}
        node.profile_variable = {}
        for pos in variable_pos:
            nuc = seq[pos]
            # divide the parent profile by the contribution coming in from this node 
            # (one could also multiply the remaining contributions, but dividing is more efficient in polytomies)
            stp = node.subtree_profile_variable.get(pos, node.subtree_profile_fixed[nuc])
            vec = node.parent.profile_variable.get(pos, node.parent.profile_fixed[nuc])/node.expQt.dot(stp)  # this is numerically tricky, need to guard against division by 0
            vec_norm = vec.sum()
            if np.isnan(vec_norm):
                import ipdb; ipdb.set_trace()
            node.outgroup_profile_variable[pos]=vec/vec_norm

            # if uncertaintly is high, keep this position
            vec = stp*node.expQt.T.dot(node.outgroup_profile_variable[pos])
            vec_norm =vec.sum()
            if vec.max()<(1-eps)*vec_norm or nuc!='ACGT'[vec.argmax()]:
                node.profile_variable[pos] = vec/vec_norm

        # report for fixed positions for each nucleotide
        node.profile_fixed = {}
        node.outgroup_profile_fixed = {}
        for nuc in 'ACGT':
            vec = node.parent.profile_fixed[nuc]/node.expQt.dot(node.subtree_profile_fixed[nuc])
            vec_norm = vec.sum()
            node.outgroup_profile_fixed[nuc] = vec/vec_norm
            node.profile_fixed[nuc] = node.subtree_profile_fixed[nuc]*node.expQt.T.dot(node.outgroup_profile_fixed[nuc])


def marginal_ancestral_sequences(tree, node, seq):
    tmp_seq = np.copy(seq)
    for pos, vec in node.profile_variable.items():
        tmp_seq[pos] = 'ACGT'[vec.argmax()]
    # just for demonstration, save to the node (could be streamed out, etc)
    node.reconstructed_seq = tmp_seq


# run the likelihood calculation
tree.logLH=0
post_order(tree, tree.root, tree.root.seq, subtree_profiles)
pre_order(tree, tree.root, tree.root.seq, outgroup_profiles)

pre_order(tree, tree.root, tree.root.seq, marginal_ancestral_sequences)

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

tt.infer_ancestral_sequences(marginal=True, reconstruct_tip_states=True)

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
print(tree.logLH, tt.sequence_LH(), tree.logLH-tt.sequence_LH())

# compare profiles between treetime and here
for tt_node, tree_node in zip(tt.tree.get_nonterminals(), tree.get_nonterminals()):
    for pos in tree_node.profile_variable:
        agree = np.abs(tree_node.profile_variable[pos] - tt_node.marginal_profile[pos]).sum()<eps2
        if not agree:
            print(pos, tree_node.profile_variable[pos], tt_node.marginal_profile[pos])

# compare reconstructed sequences between treetime and here
for tt_node, tree_node in zip(tt.tree.find_clades(), tree.find_clades()):
    if not np.all(tt.sequence(tt_node, reconstructed=True, as_string=False)==tree_node.reconstructed_seq):
        print(tt_node.name, "reconstructed_seqs don't agree")


def calculate_optimal_branch_length(tree, node, seq):
    if node == tree.root:
        return
    gtr = tree.gtr
    variable_pos = set.union(set(node.subtree_profile_variable.keys()), set(node.outgroup_profile_variable.keys()))
    t0 = node.branch_length
    t0 = len(node.muts)/len(seq)
    weights = []
    counts = []
    fixed_nuc_count = {k:v for k,v in node.nuc_composition.items()}
    for pos in variable_pos:
        nuc = seq[pos]
        tmp = np.zeros(len(gtr.eigenvals))
        for i in range(len(gtr.eigenvals)):
            tmp[i] = node.outgroup_profile_variable.get(pos, node.outgroup_profile_fixed[nuc]).dot(gtr.v_inv[i])*node.subtree_profile_variable.get(pos, node.subtree_profile_fixed[nuc]).dot(gtr.v[:,i])
        weights.append(tmp)
        counts.append(1)

        # this position is accounted for, hence we can subtract it from the count of fixed nucs 
        if nuc and nuc in 'ACGT': 
            fixed_nuc_count[nuc] -= 1
    
    for nuc in 'ACGT':
        tmp = np.zeros(len(gtr.eigenvals))
        for i in range(len(gtr.eigenvals)):
            tmp[i] = node.outgroup_profile_fixed[nuc].dot(gtr.v_inv[i])*node.subtree_profile_fixed[nuc].dot(gtr.v[:,i])
        weights.append(tmp)
        counts.append(fixed_nuc_count[nuc])

    max_iter = 10
    eps=t0*1e-3
    ii=0
    while ii<max_iter:
        coefficients = np.array([np.exp(t0*lam) for lam in gtr.eigenvals])
        f0=0
        fp0=0
        fpp0=0
        for c, w in zip(counts, weights):
            d = w.dot(coefficients)
            d1 = w.dot(gtr.eigenvals*coefficients)
            d2 = w.dot(gtr.eigenvals**2*coefficients)
            d3 = w.dot(gtr.eigenvals**3*coefficients)
            f0 += c*d1/d
            fp0 += c*(d2-d1**2/d)/d
            fpp0 += c*((d3-2*d2*d1/d+d1**3/d**2)/d-(d2-d1*d1/d)*d1/d**2)

        #t1 = t0 - 2*f0*fp0/(2*fp0**2-f0*fpp0)
        t1 = t0 - f0/fp0
        if np.abs(t1-t0)<eps:
            break
        t0 = max(0.0, t1)
        ii += 1

    node.branch_length = max(0.0, t1)
    print(node.name, node.branch_length)

tree.gtr = myGTR
for i in range(4):
    tree.logLH=0
    post_order(tree, tree.root, tree.root.seq, subtree_profiles)
    print(tree.logLH)
    pre_order(tree, tree.root, tree.root.seq, outgroup_profiles)
    pre_order(tree, tree.root, tree.root.seq, calculate_optimal_branch_length)

tree.logLH=0
post_order(tree, tree.root, tree.root.seq, subtree_profiles)
print(tree.logLH)

