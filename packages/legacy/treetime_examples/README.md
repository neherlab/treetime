## Collections of examples of phylogenetic analysis with [TreeTime](https://github.com/neherlab/treetime)

There are three different ways of using TreeTime:
 1. as a python module as part of [python scripts and programs](#treetime-in-python-scripts)
 2. [via a set of command-line tools](#command-line-usage)
 3. via the webserver at [treetime.ch](http://treetime.ch).

The documentation of TreeTime is now hosted on readthedocs: [![readthedocs](https://readthedocs.org/projects/treetime/badge/)](https://treetime.readthedocs.io/en/latest/)


![Molecular clock phylogeny of 200 NA sequences of influenza A H3N2](https://raw.githubusercontent.com/neherlab/treetime_examples/master/figures/tree_and_clock.png)

### Command line usage
The general pattern to use TreeTime to infer timetrees via the command line is
```bash
treetime --tree mytree.nwk --aln mysequences.fasta --dates metadata.csv
```
In addition, there are number of subcommands for related tasks.
Please see the following pages that explain usage for different cases:
 * [timetree estimation](https://treetime.readthedocs.io/en/latest/tutorials/timetree.html)
 * [estimation of substitution rates](https://treetime.readthedocs.io/en/latest/tutorials/clock.html)
 * [ancestral sequence reconstruction](https://treetime.readthedocs.io/en/latest/tutorials/ancestral.html)
 * [homoplasy analysis](https://treetime.readthedocs.io/en/latest/tutorials/homoplasy.html)
 * [discrete traits and mugration models](https://treetime.readthedocs.io/en/latest/tutorials/mugration.html)


### TreeTime in python scripts
A minimal script to run TreeTime would look something like this
```python
from treetime import TreeTime
from treetime.utils import parse_dates

dates=parse_dates('metadata.csv')
tt=TreeTime(tree='mytree.nwk', aln='mysequences.fasta', dates=dates)
tt.run(root='best')
```
More detailed examples can be found in the following set of scripts:
 * [rerooting and timetree inference](scripts/rerooting_and_timetrees.py)
 * [ancestral sequence reconstruction](scripts/ancestral.py)
 * [relaxed clock models](scripts/relaxed_clock.py)
 * [an analysis of 2014-2016 Ebola virus sequences](scripts/ebola.py)
