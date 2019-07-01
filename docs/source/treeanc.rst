***************************
TreeAnc class documentation
***************************

This is the core class of the TreeTime module. It stores the phylogenetic tree and implements the basic algorithms for sequence manipulation, sequence reconstruction, and branch length optimization.

The tree is stored as Bio.Phylo object. In order to facilitate the tree operations, each node of the tree is decorated with additional attributes which are set during the tree preparation. These attributes need to be updated after tree modifications. The sequences are also attached to the tree nodes. In order to save memory, the sequences are stored in the compressed form. The TreeAnc class implements methods to compress and decompress sequences.

The main purpose of the TreeAnc class is to implement standard algorithms for ancestral sequence reconstruction.
Both marginal and joint maximum likelihood reconstructions are possible.
The marginal reconstructions computes the entire distribution of the states at a given node after tracing out states at all other nodes.

The `example scripts <https://github.com/neherlab/treetime/blob/master/examples/ancestral_inference.py>`_ illustrate how to instantiate TreeAnc objects.


TreeAnc Constructor
===================

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

.. automethod:: treetime.TreeAnc.aln()

.. automethod:: treetime.TreeAnc.gtr()

.. automethod:: treetime.TreeAnc.tree()

.. automethod:: treetime.TreeAnc.leaves_lookup()


Sequence and profile manipulation
---------------------------------

.. automethod:: treetime.TreeAnc.get_mutations

.. automethod:: treetime.TreeAnc.get_reconstructed_alignment

.. automethod:: treetime.TreeAnc.make_reduced_alignment

.. automethod:: treetime.TreeAnc.expanded_sequence

.. automethod:: treetime.TreeAnc.dict_sequence


Ancestral reconstruction and tree optimization
----------------------------------------------

.. automethod:: treetime.TreeAnc.infer_ancestral_sequences

.. automethod:: treetime.TreeAnc.sequence_LH

.. automethod:: treetime.TreeAnc.optimize_tree

.. automethod:: treetime.TreeAnc.infer_gtr

.. automethod:: treetime.TreeAnc.get_reconstructed_alignment

.. automethod:: treetime.TreeAnc.get_tree_dict






