Time trees of segmented genomes with TreeKnit MCCs
--------------------------------------------------

Segments of genomes like influenza each have their own tree, but lineages share their history
in every segment unless a reassortment separated them.
`TreeKnit <https://github.com/neherlab/TreeKnit.jl>`_ infers where segment trees share their
history: the maximally compatible clades (MCCs) of each pair of trees.
``treetime arg`` uses these MCCs to infer a time tree for each segment, informed by the other
segments where they share a branch.

.. code-block:: bash

   treetime arg --trees ha_resolved.nwk na_resolved.nwk mp_resolved.nwk \
                --alignments ha.fasta na.fasta mp.fasta \
                --mccs MCCs.json --dates metadata.tsv --outdir results

``--trees`` and ``--alignments`` take two or more files, in the same order.
``--mccs`` is the ``MCCs.json`` written by TreeKnit; the pair of trees is looked up by tree
label, which TreeKnit derives from the tree file name (output trees like ``ha_resolved.nwk``
match the label ``ha``).
For two trees, the older ``MCCs.dat`` format (one MCC per line, leaves separated by commas) is
also accepted.

Each tree is analysed in turn as the focal tree. The alignments are concatenated, and every
branch of the focal tree uses the focal segment plus each other segment whose MCC contains
both ends of the branch. Ancestral reconstruction and branch lengths thus draw on all segments
that share a branch. Sequences missing from a segment are treated as missing data.

Outputs carry the index of the focal tree in ``--trees``: ``timetree_1.nexus``,
``divergence_tree_1.nexus``, ``branch_mutations_1.txt``, ``auspice_tree_1.json``, and so on.
In the Nexus files, branches with mutations are annotated with their MCC: ``mcc="i"`` with two
trees, or ``mcc_<label>="i"`` for each other tree, where ``i`` indexes that pair's MCCs in
``MCCs.json`` (``None`` outside any MCC).
