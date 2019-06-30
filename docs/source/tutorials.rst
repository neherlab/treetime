TreeTime command line usage
===========================

TreeTime implements a command line interface (for details, see `Command-line API <commands.html>`_) that allow estimation of time scaled phylogenies, ancestral reconstruction, and analysis of temporal signal in alignments.
The command interface is organized as the main command performing time-tree estimation as

.. code:: bash

  treetime --tree tree_file --aln alignment --dates dates.tsv

with other functionalities available as subcommands

.. code:: bash

  treetime {ancestral, clock, homoplasy, mugration, version}


TreeTime can use full alignments in `fasta` or `phylip` format or work of VCF files.
For each of the different subcommands, we prepared tutorials listed below.
These tutorials use example data provided in the github-repository
`github.com/neherlab/treetime_examples <https://github.com/neherlab/treetime_examples>`_.


.. toctree::
    :maxdepth: 1

    tutorials/timetree
    tutorials/ancestral
    tutorials/clock
    tutorials/mugration
    tutorials/homoplasy