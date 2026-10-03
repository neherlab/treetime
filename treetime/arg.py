import itertools
import json
import os

import numpy as np


def read_mccs(mcc_file, tree_labels=None):
    """Read the maximally compatible clades (MCCs) inferred by TreeKnit.

    Two formats are accepted:
    - JSON (``MCCs.json``, TreeKnit >= 0.5):
      ``{"MCC_dict": {"1": {"trees": [label_a, label_b], "mccs": [[leaf, ...], ...]}, ...}}``
      with one entry per pair of trees.
    - text (``MCCs.dat``): one MCC per line, leaves separated by commas.

    Args:
        mcc_file (str): file name
        tree_labels (list, optional): labels of the two trees, used to select the pair from a
            JSON file with several pairs. TreeKnit labels trees by their file name without
            extension. Ignored for text files and for JSON files with a single pair.

    Returns:
        list: MCCs as lists of leaf names
    """
    with open(mcc_file) as fh:
        content = fh.read()

    if not content.lstrip().startswith('{'):
        return [line.strip().split(',') for line in content.splitlines() if line.strip()]

    pairs = list(json.loads(content)['MCC_dict'].values())
    if not pairs:
        raise ValueError(f'no tree pairs in {mcc_file}')
    if tree_labels is not None:
        for pair in pairs:
            if set(pair['trees']) == set(tree_labels):
                return pair['mccs']
    if len(pairs) == 1:
        return pairs[0]['mccs']
    available = ', '.join('/'.join(p['trees']) for p in pairs)
    raise ValueError(
        f'{mcc_file} contains several tree pairs ({available}); none matches the trees '
        f'{"/".join(tree_labels or [])}. TreeKnit labels trees by file name without extension.'
    )


def tree_label(tree_file):
    """Label TreeKnit gives the tree in `tree_file`: the file name without extension, and without
    the suffixes of TreeKnit's output trees (`_resolved`, `_liberal_resolved`, `_imputed`), so that
    trees written by TreeKnit match the labels in its MCC file."""
    label = os.path.splitext(os.path.basename(tree_file))[0]
    for suffix in ('_liberal_resolved', '_resolved', '_imputed'):
        if label.endswith(suffix):
            return label[: -len(suffix)]
    return label


def parse_arg(trees, alignments, mcc_file, fill_overhangs=True):
    """Read the input of `treetime arg`: K >= 2 trees, one alignment per tree, and the MCCs
    inferred by TreeKnit.

    Args:
        trees (list): file names of the trees
        alignments (list): file names of the alignments, in the same order as `trees`
        mcc_file (str): MCC file written by TreeKnit (JSON or text, see `read_mccs`). Text files
            contain a single MCC set and can only be used with two trees.
        fill_overhangs (bool, optional): fill terminal gaps of alignments before concatenating.
            Defaults to True.

    Returns:
        dict with
        - 'trees': the trees (Bio.Phylo)
        - 'labels': tree labels as used by TreeKnit (see `tree_label`)
        - 'alignment': concatenation of all alignments over all leaves of all trees; sequences
          missing from a segment are filled with N (uninformative)
        - 'masks': one array per tree, 1 at the positions of its segment in the concatenation
        - 'combined_mask': array of ones over the whole concatenation
        - 'MCCs': {(i, j): MCCs of trees i and j} for all pairs i < j
    """
    from Bio import AlignIO, Phylo, Seq
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord

    from treetime.seq_utils import seq2array

    if len(trees) < 2 or len(trees) != len(alignments):
        raise ValueError('treetime arg needs at least two trees and one alignment per tree')
    labels = [tree_label(t) for t in trees]
    if len(set(labels)) != len(labels):
        raise ValueError(f'tree labels must be unique, got {labels}')

    parsed_trees = [Phylo.read(t, 'newick') for t in trees]
    all_leaves = sorted(set.union(*[{x.name for x in t.get_terminals()} for t in parsed_trees]))

    # read MCCs of every pair of trees as lists of taxon names
    MCCs = {}
    for i, j in itertools.combinations(range(len(trees)), 2):
        MCCs[i, j] = read_mccs(mcc_file, tree_labels=[labels[i], labels[j]])

    # read alignments with edge-modified sequences; missing sequences are all N
    segments = []
    for aln_file in alignments:
        aln = {s.id: ''.join(seq2array(s, fill_overhangs=fill_overhangs)) for s in AlignIO.read(aln_file, 'fasta')}
        length = len(next(iter(aln.values())))
        missing = [leaf for leaf in all_leaves if leaf not in aln]
        if missing:
            print(f'WARNING: {len(missing)} leaves without a sequence in {aln_file} (treated as missing data)')
        segments.append((aln, length))

    aln_combined = MultipleSeqAlignment(
        [
            SeqRecord(Seq.Seq(''.join(aln.get(leaf, 'N' * length) for aln, length in segments)), id=leaf)
            for leaf in all_leaves
        ]
    )

    # masks: the positions of each segment in the concatenation
    total = sum(length for _, length in segments)
    masks = []
    start = 0
    for _, length in segments:
        mask = np.zeros(total)
        mask[start : start + length] = 1
        masks.append(mask)
        start += length

    return {
        'trees': parsed_trees,
        'labels': labels,
        'alignment': aln_combined,
        'masks': masks,
        'combined_mask': np.ones(total),
        'MCCs': MCCs,
    }


def setup_arg(
    T,
    aln,
    focal_mask,
    other_masks,
    other_mccs,
    dates,
    gtr='JC69',
    verbose=0,
    fill_overhangs=True,
    reroot=True,
    fixed_clock_rate=None,
    alphabet='nuc',
    **kwargs,
):
    """Construct a TreeTime object for the tree of a focal segment, informed by other segments.

    Each branch of the focal tree is annotated with its MCC for every other segment
    (`n.mccs[label]`) and gets a mask over the concatenated alignment: the focal segment, plus
    every other segment that shares the branch, i.e. whose MCC contains both ends of the branch.
    Ancestral reconstruction and branch length optimization then use the other segments only
    where they share the focal segment's history.

    Args:
        T (str, Bio.Phylo.Tree): tree of the focal segment
        aln (Bio.Align.MultipleSeqAlignment): concatenated multiple sequence alignment
        focal_mask (np.array): 1 at the positions of the focal segment
        other_masks (dict): label of another segment -> mask of its positions
        other_mccs (dict): label of another segment -> MCCs of the focal and that segment
        dates (dict): sampling dates
        gtr (str, optional): GTR model. Defaults to 'JC69'.
        verbose (int, optional): verbosity. Defaults to 0.
        fill_overhangs (bool, optional): treat terminal gap as missing. Defaults to True.
        reroot (bool, optional): reroot the tree. Defaults to True.

    Returns:
        TreeTime: TreeTime instance
    """
    from treetime import TreeTime

    tt = TreeTime(
        dates=dates,
        tree=T,
        aln=aln,
        gtr=gtr,
        alphabet=alphabet,
        verbose=verbose,
        fill_overhangs=fill_overhangs,
        keep_node_order=True,
        compress=False,
        **kwargs,
    )

    if reroot:
        tt.reroot('least-squares', force_positive=True, clock_rate=fixed_clock_rate)

    enforce_min_branch_length(tt.tree, tt.one_mutation)
    annotate_mccs(tt.tree, other_mccs)
    assign_masks(tt.tree, focal_mask, other_masks)
    return tt


def annotate_mccs(tree, mccs_by_label):
    """Set `n.mccs = {label: MCC index or None}` on every node of `tree`, one entry per other
    segment, from that segment's MCCs with the focal segment. With a single other segment,
    `n.mcc` is set as well.
    """
    for n in tree.find_clades():
        n.mccs = {}
    for label, mccs in mccs_by_label.items():
        leaf_to_mcc = {leaf: mi for mi, mcc in enumerate(mccs) for leaf in mcc}
        for n, m in map_mccs(tree, leaf_to_mcc).items():
            n.mccs[label] = m
    for n in tree.find_clades():
        n.mcc = next(iter(n.mccs.values())) if len(n.mccs) == 1 else None


def assign_masks(tree, focal_mask, other_masks):
    """Set `n.mask` on every node: the focal segment, plus each other segment whose MCC contains
    both the node and its parent (the branch is shared). Requires `annotate_mccs`."""
    for n in tree.find_clades():
        mask = focal_mask.copy()
        if n.up is not None:
            for label, segment_mask in other_masks.items():
                m = n.mccs.get(label)
                if m is not None and n.up.mccs.get(label) == m:
                    mask = mask + segment_mask
        n.mask = mask


def enforce_min_branch_length(tree, one_mutation=1e-4):
    """Give every branch at least half a mutation of length."""
    for n in tree.find_clades():
        if n.branch_length is not None:
            n.branch_length = max(0.5 * one_mutation, n.branch_length)


def map_mccs(tree, mcc_map):
    """MCC of every node of `tree` (Fitch-like two-pass assignment).

    Args:
        tree (Bio.Phylo.Tree): tree
        mcc_map (dict): leaf name -> MCC index. Leaves without an entry (e.g. absent from the
            other segment) constrain nothing.

    Returns:
        dict: node -> MCC index, or None for nodes not inside an MCC
    """
    # upward pass: sets of MCCs consistent with the subtree; None means unconstrained
    up = {}
    for n in tree.find_clades(order='postorder'):
        if n.is_terminal():
            up[n] = {mcc_map[n.name]} if n.name in mcc_map else None
            continue
        sets = [up[c] for c in n if up[c] is not None]
        if not sets:
            up[n] = None
            continue
        common = set.intersection(*sets)
        up[n] = common if (common or n == tree.root) else set.union(*sets)

    # downward pass
    mcc = {}
    for n in tree.find_clades(order='preorder'):
        s = up[n]
        if s is None:
            mcc[n] = None
        elif n == tree.root:  # children share an MCC
            mcc[n] = next(iter(s)) if s else None
        elif mcc[n.up] in s:  # parent MCC continues into this node
            mcc[n] = mcc[n.up]
        elif len(s) == 1:  # root of an MCC
            mcc[n] = next(iter(s))
        else:  # not inside an MCC
            mcc[n] = None
    return mcc


def assign_mccs(tree, mcc_map, one_mutation=1e-4):
    """Assign MCCs to all terminal and internal branches of the tree (`n.mcc`) and give every
    branch at least half a mutation of length.

    Args:
        tree (Bio.Phylo.Tree): tree
        mcc_map (dict): map from leaf to mcc
        one_mutation (float, optional): minimal length of branches. Defaults to 1e-4.
    """
    enforce_min_branch_length(tree, one_mutation)
    for n, m in map_mccs(tree, mcc_map).items():
        n.mcc = m
