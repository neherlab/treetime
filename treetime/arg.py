from matplotlib.pyplot import fill
import numpy as np
import json
import itertools

def get_tree_names(tree_nwk_files):
    tree_names = []
    for file in tree_nwk_files:
        file_name = file.split("/")[-1].split(".")[0]
        file_name = file_name.replace("_resolved", "").replace("resolved", "")
        tree_names.append(file_name)
    if len(set(tree_names)) != len(tree_nwk_files):
        #tree name identifiers are not unique
        raise Exception("Error: Tree names must be unique, see TreeKnit output format.")
    return tree_names

def get_MCC_dict(MCC_file):   
    f = open(MCC_file)
    data = json.load(f)
    MCC_dict = {}
    for key in data["MCC_dict"]:
        MCC_dict[frozenset(data["MCC_dict"][key]["trees"])] = data["MCC_dict"][key]["mccs"]

    return MCC_dict

def get_mask_dict(length_segments, tree_names):
    pos_list = [0]
    for l in length_segments:
        new = pos_list[-1] + l
        pos_list.append(new)
    mask = {}
    no_trees = len(tree_names)
    for r in range(1,(no_trees+1)):
        combos = itertools.combinations(range(1, (no_trees+1)), r)
        for comb in combos:
            new_mask = np.zeros(sum(length_segments))
            for c in comb:
                new_mask[pos_list[c-1]:pos_list[c]] = 1
            mask[frozenset([tree_names[c-1] for c in comb])] = new_mask
    return mask

def parse_arg(tree_files, aln_files, MCC_file, fill_overhangs=True):
#def parse_arg(tree_file_dir, aln_file_dir, MCC_file, tree_names, fill_overhangs=True):
    """parse the output of TreeKnit and return a file structure to be
    further consumed by TreeTime

    Args:
        tree_files (str): file names of trees
        aln_files (str): file names of alignments MUST be in the same order as tree_files
        MCC_file (str): name of mcc file
        fill_overhangs (bool, optional): fill terminal gaps of alignmens before concatenating. Defaults to True.

    Returns:
        dict: dictionary containing the two trees, the concatenated alignment, full and segment masks, and the MCCs
    """
    from Bio import Phylo, AlignIO, Seq
    from Bio.Align import MultipleSeqAlignment
    from treetime.seq_utils import seq2array

    # read MCCs as lists of taxon names
    MCC_dict = get_MCC_dict(MCC_file)

    # read trees, assert trees are named consistently
    tree_names = get_tree_names(tree_files)
    trees_in_dict = set().union(*MCC_dict.keys())
    assert(all([k in trees_in_dict for k in tree_names]))
    trees_dict = {}
    for i in range(0, len(tree_files)):
        trees_dict[tree_names[i]] = Phylo.read(tree_files[i], 'newick')

    # determine common terminal nodes
    all_leaves = set.intersection(*[set([x.name for x in t.get_terminals()]) for (k, t) in trees_dict.items()])

    # read alignments and construct edge modified sequence arrays
    alignments = [{s.id:s for s in AlignIO.read(aln, 'fasta')} for aln in aln_files]
    for aln in alignments:
        for s,seq in aln.items():
            seqstr = "".join(seq2array(seq, fill_overhangs=fill_overhangs))
            seq.seq = Seq.Seq(seqstr)

    # construct concatenated alignment
    aln_combined = []
    for leaf in all_leaves:
        concat_seq = alignments[1][leaf]
        for a in range(1, len(alignments)):
            concat_seq += alignments[a][leaf]
        seq = concat_seq
        seq.id = leaf
        aln_combined.append(seq)

    # construct masks for the concatenation and the two segments
    l = [len(a[leaf]) for a in alignments]
    masks = get_mask_dict(l, tree_names)

    return {"MCCs_dict": MCC_dict, "trees_dict":trees_dict, "alignment":MultipleSeqAlignment(aln_combined),
            "masks_dict":masks}

def setup_arg(trees_dict, alignments, dates, MCCs_dict, masks_dict, tree_name, gtr='JC69',
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

    T= trees_dict[tree_name] ##desired tree

    ##get list of MCCs of all trees with T and the order of these trees
    MCCs = []
    tree_order = {}
    i = 0
    for t in trees_dict.keys():
        if t != tree_name:
            tree_order[i] = t
            MCCs.append(MCCs_dict[frozenset([tree_name, t])])
            i +=1

    tt = TreeTime(dates=dates, tree=T,
            aln=alignments, gtr=gtr, alphabet=alphabet, verbose=verbose,
            fill_overhangs=fill_overhangs, keep_node_order=True,
            compress=False, **kwargs)


    if reroot:
        tt.reroot("least-squares", force_positive=True, clock_rate=fixed_clock_rate)

    # make a lookup for the MCCs and assign to tree
    leaf_to_MCC = get_mcc_map(MCCs)

    assign_all_mccs(tt.tree, len(MCCs), leaf_to_MCC, tt.one_mutation)

    # assign masks to branches whenever child and parent are in the same MCC
    for n in tt.tree.find_clades():
        shared = [(n.mcc[pos] is not None) and n.up and n.up.mcc[pos]==n.mcc[pos] for pos in range(len(MCCs))]
        ##use tree_order to convert position in MCC list to tree_names and see which trees share this branch and assign a proper mask
        pos_shared = [tree_order[i] for i, x in enumerate(shared) if x]
        pos_shared.append(tree_name)
        n.mask = masks_dict[frozenset(pos_shared)]


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

def get_mcc_map(MCCs_list):
    # make a lookup for the MCCs and assign to trees
    leaf_to_MCC = {}
    for MCCs in MCCs_list:
        for mi,mcc in enumerate(MCCs):
            for leaf in mcc:
                if leaf not in leaf_to_MCC:
                    leaf_to_MCC[leaf] = [mi]
                else:
                    leaf_to_MCC[leaf].append(mi)
    return leaf_to_MCC

def assign_all_mccs(tree, len_tree_list, mcc_map, one_mutation=1e-4):
    for leaf in tree.get_terminals():
        leaf.child_mccs = [set([mcc_map[leaf.name][pos]]) for pos in range(len_tree_list)]
        leaf.mcc = mcc_map[leaf.name]
        leaf.branch_length = max(0.5*one_mutation, leaf.branch_length)
    # reconstruct MCCs with Fitch algorithm
    for n in tree.get_nonterminals(order='postorder'):
        common_mccs = [set.intersection(*[c.child_mccs[pos] for c in n]) for pos in range(len_tree_list)]
        n.branch_length = max(0.5*one_mutation, n.branch_length)
        n.child_mccs = []
        for pos in range(len_tree_list):
            if len(common_mccs[pos]):
                n.child_mccs.append(common_mccs[pos])
            else:
                n.child_mccs.append(set.union(*[c.child_mccs[pos] for c in n]))
    mcc_intersection = [set.intersection(*[c.child_mccs[pos] for c in tree.root]) for pos in range(len_tree_list)]
    tree.root.mcc = []
    for pos in range(len_tree_list):
        if len(mcc_intersection[pos]):
            tree.root.mcc.append(list(mcc_intersection[pos])[0])
        else:
            tree.root.mcc.append(None)
    for n in tree.get_nonterminals(order='preorder'):
        if n==tree.root:
            continue
        else:
            n.mcc = []
            for pos in range(len_tree_list):
                if n.up.mcc[pos] in n.child_mccs[pos]: # parent MCC part of children -> that is the MCC
                    n.mcc.append(n.up.mcc[pos])
                elif len(n.child_mccs[pos])==1:  # child is an MCC
                    n.mcc.append(list(n.child_mccs[pos])[0])
                else: # no unique child MCC and no match with parent -> not part of an MCCs
                    n.mcc.append(None)