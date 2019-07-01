.. TreeTime documentation master file, created by
   sphinx-quickstart on Mon Jul 31 11:44:07 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TreeTime: time-tree and ancestral sequence inference
====================================================

.. image:: https://travis-ci.org/neherlab/treetime.svg?branch=master
   :target: https://travis-ci.org/neherlab/treetime

.. image:: https://anaconda.org/bioconda/treetime/badges/installer/conda.svg
   :target: https://anaconda.org/bioconda/treetime

TreeTime provides routines for ancestral sequence reconstruction and inference of molecular-clock phylogenies, i.e., a tree where all branches are scaled such that the positions of terminal nodes correspond to their sampling times and internal nodes are placed at the most likely time of divergence.

To optimize the likelihood of time-scaled phylogenies, TreeTime uses an iterative approach that first optimizes branch lengths of the tree given the sequence data and date constraints, and then optimizes coalescent tree priors, relaxed clock parameters, or resolves polytomies.
This cycle is repeated a few times.
The only topology optimization are (optional) resolution of polytomies in a way that is most (approximately) consistent with the sampling time constraints on the tree.

The code is hosted on `github.com/neherlab/treetime <https://github.com/neherlab/treetime>`_.

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   tutorials
   commands
   APIdoc


.. image:: https://raw.githubusercontent.com/neherlab/treetime_examples/master/figures/tree_and_clock.png

Features
--------

  * ancestral sequence reconstruction (marginal and joint maximum likelihood)

  * molecular clock tree inference (marginal and joint maximum likelihood)

  * inference of GTR models

  * rerooting to maximize temporal signal and optimize the root-to-tip distance vs time relationship

  * simple phylodynamic analysis such as coalescent model fits


Developer info
--------------
  - Source code on github at https://github.com/neherlab/treetime

  - Copyright and License: Pavel Sagulenko, Emma Hodcroft, and Richard Neher, MIT Licence

  - References

    * `TreeTime: Maximum-likelihood phylodynamic analysis <https://academic.oup.com/ve/article/4/1/vex042/4794731>`_ by Pavel Sagulenko, Vadim Puller and Richard A Neher. Virus Evolution.

    * `NextStrain: real-time tracking of pathogen evolution <https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty407/5001388>`_ by James Hadfield et al. Bioinformatics.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

