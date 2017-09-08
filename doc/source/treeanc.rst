***************************
TreeAnc class documentation
***************************

This is the basic class of the TreeTime module. Its purpose to store the phylogenetic tree, and to implement the basic algorithms for sequence manipulation, reconstruction, and branch length optimization.

The tree is stored as Bio.Phylo object. In order to facilitate the tree operations, each node of the tree gets some additional attributes, which set during the tree preparation. These attributes should be updated properly after tree modifications. The sequences are also attached to the tree nodes. In order to save space, the sequences are stored in the compressed form, so the TreeAnc class has methods to compress and decompress sequences.

On top of that, the TreeAnc class realizes some standard popular algorithms for ancestral sequence reconstruction, maximum-likelihood algorithm to optimize a branch length knowing sequences on each side of the tree branch.

There are some examples given to show how to instantiate objects of TreeAnc type, and how to run ancestral sequence reconstruction and branch lengths optimization.


TreeAnc Docstring
=================

.. autoclass:: treetime.TreeAnc
    :members: __init__

TreeAnc methods
===============

Basic functions, utilities, properties
--------------------------------------

.. automethod:: treetime.TreeAnc.prepare_tree

.. automethod:: treetime.TreeAnc.prune_short_branches

.. automethod:: treetime.TreeAnc.set_gtr

.. automethod:: treetime.TreeAnc.logger

.. automethod:: treetime.TreeAnc.aln

.. automethod:: treetime.TreeAnc.gtr

.. automethod:: treetime.TreeAnc.tree

.. automethod:: treetime.TreeAnc.leaves_lookup


Sequence and profile manipulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: treetime.TreeAnc.get_mutations

.. automethod:: treetime.TreeAnc.get_reconstructed_alignment

.. automethod:: treetime.TreeAnc.make_reduced_alignment

.. automethod:: treetime.TreeAnc.expanded_sequence


Ancestral reconstruction and tree optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: treetime.TreeAnc.reconstruct_anc

.. automethod:: treetime.TreeAnc.infer_ancestral_sequences

.. automethod:: treetime.TreeAnc.ancestral_likelihood

.. automethod:: treetime.TreeAnc.optimal_branch_length

.. automethod:: treetime.TreeAnc.optimize_branch_length

.. automethod:: treetime.TreeAnc.optimize_seq_and_branch_len

.. automethod:: treetime.TreeAnc.optimize_sequences_and_branch_length

.. automethod:: treetime.TreeAnc.infer_gtr








