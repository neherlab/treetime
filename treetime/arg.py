from matplotlib.pyplot import fill
import numpy as np

def parse_arg(tree1, tree2, aln1, aln2, MCC_file, fill_overhangs=True):
    """parse the output of TreeKnit and return a file structure to be
    further consumed by TreeTime

    Args:
        tree1 (str): file name of tree1
        tree2 (str): file name of tree2
        aln1 (str): file name of alignment 1
        aln2 (str): file name of alignment 2
        MCC_file (str): name of mcc file
        fill_overhangs (bool, optional): fill terminal gaps of alignmens before concatenating. Defaults to True.

    Returns:
        dict: dictionary containing the two trees, the concatenated alignment, full and segment masks, and the MCCs
    """
    from Bio import Phylo, AlignIO, Seq
    from Bio.Align import MultipleSeqAlignment
    from treetime.seq_utils import seq2array

    # read trees and determine common terminal nodes
    t1 = Phylo.read(tree1, 'newick')
    t2 = Phylo.read(tree2, 'newick')
    all_leaves = set.intersection(set([x.name for x in t1.get_terminals()]), set([x.name for x in t2.get_terminals()]))

    # read MCCs as lists of taxon names
    MCCs = []
    with open(MCC_file) as fh:
        for line in fh:
            if line.strip():
                MCCs.append(line.strip().split(','))

    # read alignments and construct edge modified sequence arrays
    a1 = {s.id:s for s in AlignIO.read(aln1, 'fasta')}
    a2 = {s.id:s for s in AlignIO.read(aln2, 'fasta')}
    for aln in [a1,a2]:
        for s,seq in aln.items():
            seqstr = "".join(seq2array(seq, fill_overhangs=fill_overhangs))
            seq.seq = Seq.Seq(seqstr)

    # construct concatenated alignment
    aln_combined = []
    for leaf in all_leaves:
        seq = a1[leaf] + a2[leaf]
        seq.id = leaf
        aln_combined.append(seq)

    # construct masks for the concatenation and the two segments
    l1 = len(a1[leaf])
    l2 = len(a2[leaf])
    combined_mask = np.ones(l1 + l2)
    mask1 = np.zeros(l1 + l2)
    mask2 = np.zeros(l1 + l2)
    mask1[:l1] = 1
    mask2[l1:] = 1

    return {"MCCs": MCCs, "trees":[t1,t2], "alignment":MultipleSeqAlignment(aln_combined),
            "masks":[mask1,mask2], "combined_mask":combined_mask}

def setup_arg(T, aln, total_mask, segment_mask, dates, MCCs, gtr='JC69',
            verbose=0, fill_overhangs=True, reroot=True, fixed_clock_rate=None, alphabet='nuc', **kwargs):
    """construct a TreeTime object with the appropriate masks on each node
    for branch length optimization with full or segment only alignment.

    Args:
        T (str, Bio.Phylo.Tree): tree of focal segment
        aln (Bio.Align.MultipleSeqAlignment): Concatenated multiple sequence alignment
        total_mask (np.array): boolean array that is true for the entire sequence
        segment_mask (np.array): boolean array that is true only for the focal segment
        dates (dict): sampling dates
        MCCs (list): list of MCCs
        gtr (str, optional): GTR model. Defaults to 'JC69'.
        verbose (int, optional): verbosity. Defaults to 0.
        fill_overhangs (bool, optional): treat terminal gap as missing. Defaults to True.
        reroot (bool, optional): reroot the tree. Defaults to True.

    Returns:
        TreeTime: TreeTime instance
    """
    from treetime import TreeTime

    tt = TreeTime(dates=dates, tree=T,
            aln=aln, gtr=gtr, alphabet=alphabet, verbose=verbose,
            fill_overhangs=fill_overhangs, keep_node_order=True,
            compress=False, **kwargs)


    if reroot:
        tt.reroot("least-squares", force_positive=True, clock_rate=fixed_clock_rate)

    # make a lookup for the MCCs and assign to tree
    leaf_to_MCC = {}
    for mi,mcc in enumerate(MCCs):
        for leaf in mcc:
            leaf_to_MCC[leaf] = mi

    assign_mccs(tt.tree, leaf_to_MCC, tt.one_mutation)

    # assign masks to branches whenever child and parent are in the same MCC
    for n in tt.tree.find_clades():
        if (n.mcc is not None) and n.up and n.up.mcc==n.mcc:
            n.mask = total_mask
        else:
            n.mask = segment_mask

    return tt


def assign_mccs(tree, mcc_map, one_mutation=1e-4):
    """Assign MCCs to all terminal and internal branches of the tree.

    Args:
        tree (Bio.Phylo.Tree): tree
        mcc_map (dict): map from leaf to mcc
        one_mutation (float, optional): minimal length of branches. Defaults to 1e-4.
    """
    # assign MCCs to leaves
    for leaf in tree.get_terminals():
        leaf.child_mccs = set([mcc_map[leaf.name]])
        leaf.mcc = mcc_map[leaf.name]
        leaf.branch_length = max(0.5*one_mutation, leaf.branch_length)

    # reconstruct MCCs with Fitch algorithm
    for n in tree.get_nonterminals(order='postorder'):
        common_mccs = set.intersection(*[c.child_mccs for c in n])
        n.branch_length = max(0.5*one_mutation, n.branch_length)
        if len(common_mccs):
            n.child_mccs = common_mccs
        else:
            n.child_mccs = set.union(*[c.child_mccs for c in n])

    mcc_intersection = set.intersection(*[c.child_mccs for c in tree.root])
    if len(mcc_intersection):
        tree.root.mcc = list(mcc_intersection)[0]
    else:
        tree.root.mcc = None

    for n in tree.get_nonterminals(order='preorder'):
        if n==tree.root:
            continue
        else:
            if n.up.mcc in n.child_mccs: # parent MCC part of children -> that is the MCC
                n.mcc = n.up.mcc
            elif len(n.child_mccs)==1:  # child is an MCC
                n.mcc = list(n.child_mccs)[0]
            else: # no unique child MCC and no match with parent -> not part of an MCCs
                n.mcc = None
