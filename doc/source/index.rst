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
TreeTime is designed to be part of python workflows, but we have exposed a number of standard
tasks via a command-line interface.
The TreeTime command-line tool is called :code:`treetime`.
Examples and documentation of the command-line interface can be found in the github repo https://github.com/neherlab/treetime_examples.
In its standard mode, it will take a tree, an alignment, and file with dates as input and estimate a time-scaled phylogeny.
The full set of options are available via :code:`treetime -h`.


Subcommand :code:`treetime ancestral`
-------------------------------------
This subcommand reconstructs ancestral sequences and maps mutations to the tree.
It produces an alignment file containing inferred ancestral sequences and a tree file
with mutations included as comments. The inferred GTR model is written to stdout.

Subcommand :code:`treetime homoplasy`
-------------------------------------
Reconstructs ancestral sequences and maps mutations to the tree.
The tree is then scanned for homoplasies. An excess number of homoplasies
might suggest contamination, recombination, culture adaptation or similar.
Results are printed to stdout.


Subcommand :code:`treetime clock`
---------------
Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree.
It will reroot the tree to maximize the clock-like
signal and recalculate branch length unless run with :code:`--keep_root`.

Subcommand :code:`treetime mugration`
-------------------------------------
Reconstructs discrete ancestral states, for example geographic location, host, or similar.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

