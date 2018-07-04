.. TreeTime documentation master file, created by
   sphinx-quickstart on Mon Jul 31 11:44:07 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TreeTime documentation
====================================

TreeTime is organized as a hierarchy of classes. The GTR class implements sequence evolution models, TreeAnc does ancestral sequence reconstruction, ClockTree implements time tree inference for a fixed tree topology, while TreeTime provides convenient wrapper functions and additional functionality to manipulate the tree (e.g. rerooting and polytomy resolution).

.. toctree::
   :maxdepth: 2
   :hidden:

   gtr
   treeanc
   clock_tree
   treetime
   vcf_utils
   seq_utils


.. automodule:: treetime


:doc:`GTR class<gtr>`
---------------------

:doc:`TreeAnc class<treeanc>`
------------------------------

:doc:`ClockTree class<clock_tree>`
----------------------------------

:doc:`TreeTime class<treetime>`
-------------------------------


Utility code
============

:doc:`VCF tools<vcf_utils>`
-------------------------------

:doc:`Seq tools<seq_utils>`
-------------------------------


Command-line functions
======================
TreeTime is designed to be part of python workflows, but for a number of standard
tasks we have implemented scripts that can be called from the command line like
regular linux programs.

homoplasy_scanner
-----------------
Reconstructs ancestral sequences and maps mutations to the tree.
The tree is then scanned for homoplasies. An excess number of homoplasies
might suggest contamination, recombination, culture adaptation or similar.
Results are printed to stdout.

ancestral_reconstruction
------------------------
Reconstructs ancestral sequences and maps mutations to the tree.
Produces an alignment file containing inferred ancestral sequences and a tree file
with mutations included as comments. The inferred GTR model is written to stdout.

temporal_signal
---------------
Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree.
It will reroot the tree to maximize the clock-like
signal and recalculate branch length unless run with --keep_root.

timetree_inference
------------------
Reconstructs ancestral sequences and infers a molecular clock tree. Produces
an alignment file containing inferred ancestral sequences and a tree file with
mutations included as comments, and prints the molecular clock and inferred
GTR model.

mugration
---------
Reconstructs discrete ancestral states, for example geographic location, host, or similar.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

