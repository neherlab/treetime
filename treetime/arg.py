from matplotlib.pyplot import fill
import numpy as np

def parse_arg(tree1, tree2, aln1, aln2, MCC_file, fill_overhangs=True):
    from Bio import Phylo, AlignIO, Seq
    from Bio.Align import MultipleSeqAlignment
    from treetime.seq_utils import seq2array

    t1 = Phylo.read(tree1, 'newick')
    t2 = Phylo.read(tree2, 'newick')

    MCCs = []
    with open(MCC_file) as fh:
        for line in fh:
            if line.strip():
                MCCs.append(line.strip().split(','))

    a1 = {s.id:s for s in AlignIO.read(aln1, 'fasta')}
    a2 = {s.id:s for s in AlignIO.read(aln2, 'fasta')}
    all_leaves = set.intersection(set([x.name for x in t1.get_terminals()]), set([x.name for x in t2.get_terminals()]))
    for aln in [a1,a2]:
        for s,seq in aln.items():
            seqstr = "".join(seq2array(seq, fill_overhangs=fill_overhangs))
            seq.seq = Seq.Seq(seqstr)

    aln_combined = []
    for leaf in all_leaves:
        seq = a1[leaf] + a2[leaf]
        seq.id = leaf
        aln_combined.append(seq)

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
              verbose=0, fill_overhangs=True, reroot=True):
    from treetime import TreeTime
    from collections import defaultdict

    tt = TreeTime(dates=dates, tree=T,
                  aln=aln, gtr=gtr, verbose=verbose,
                  fill_overhangs=fill_overhangs, keep_node_order=True,
                  compress=False)


    if reroot:
        tt.reroot("least-squares", force_positive=True)

    leaf_to_MCC = {}
    for mi,mcc in enumerate(MCCs):
        for leaf in mcc:
            leaf_to_MCC[leaf] = mi

    for leaf in tt.tree.get_terminals():
        leaf.child_mccs = set([leaf_to_MCC[leaf.name]])
        leaf.mcc = leaf_to_MCC[leaf.name]
        leaf.branch_length = max(0.5*tt.one_mutation, leaf.branch_length)

    for n in tt.tree.get_nonterminals(order='postorder'):
        common_mccs = set.intersection(*[c.child_mccs for c in n])
        n.branch_length = max(0.5*tt.one_mutation, n.branch_length)
        if len(common_mccs):
            n.child_mccs = common_mccs
        else:
            n.child_mccs = set.union(*[c.child_mccs for c in n])

    mcc_intersection = set.intersection(*[c.child_mccs for c in tt.tree.root])
    if len(mcc_intersection):
        tt.tree.root.mcc = list(mcc_intersection)[0]
    else:
        tt.tree.root.mcc = None

    for n in tt.tree.get_nonterminals(order='preorder'):
        if n==tt.tree.root:
            continue
        else:
            if n.up.mcc in n.child_mccs:
                n.mcc = n.up.mcc
            elif len(n.child_mccs)==1:
                n.mcc = list(n.child_mccs)[0]
            else:
                n.mcc = None

    for n in tt.tree.find_clades():
        if (n.mcc is not None) and n.up and n.up.mcc==n.mcc:
            n.mask = total_mask
        else:
            n.mask = segment_mask

    return tt
