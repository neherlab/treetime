[![Build Status](https://travis-ci.org/neherlab/treetime.svg?branch=master)](https://travis-ci.org/neherlab/treetime)

# TreeTime: Maximum-Likelihood dating and ancestral inference for phylogenetic trees

## Overview

TreeTime provides routines for ancestral sequence reconstruction and the inference of molecular-clock phylogenies, i.e., a tree where all branches are scaled such that the locations of terminal nodes correspond to their sampling times and internal nodes are placed at the most likely time of divergence.

TreeTime aims at being a compromise between sophisticated probabilistic models of evolution and fast heuristics. It implements GTR models of ancestral inference and branch length optimization, but takes the tree topology as given.
The only topology optimization are resolution of polytomies in a way that is most (approximately) consistent with the sampling time constraints on the tree.
The package is designed to be used as a stand-alone tool as well as a module plugged in a bigger phylogenetic tool.

#### Features
* ancestral sequence reconstruction (marginal and joint maximum likelihood)
* molecular clock tree inference (marginal and joint maximum likelihood)
* inference of GTR models
* rerooting to obtain best root-to-tip regression
* auto-correlated relaxed molecular clock (with normal prior)


## Getting started

### Installation and prerequisites

* The package depends on several python libraries:
    - numpy, SciPy: for all kind of mathematical operations as matrix operations, numerical integration, interpolation, minimization, etc.

    - BioPython: for parsing multiple sequence alignments and all phylogenetic functionality
  If you do not have these libraries, you can install them by typing in the terminal:
    ```bash
    $pip install numpy scipy biopython
    ```

* To install the package, run `setup.py` script from the terminal:
    ```bash
    $python setup.py install
    ```

You might need root privileges for system wide installation. Alternatively, you can simply use it TreeTime locally without installation. In this case, just download and unpack it, and then add the TreeTime folder to your $PYTHONPATH.


### Basic usage

* Ancestral sequence reconstruction:

    ```python
    from treetime import TreeAnc
    ta = TreeAnc(tree='my_tree.nwk', aln='my_seqs.nwk', gtr='JC69')
    ta.reconstruct_anc(method = 'ml', infer_gtr=True, marginal=False)
    ```
  Every node of `ta.tree` now has a `node.sequence` attached. With the Optimal arguments to `infer_gtr=True`, a maximum likelihood GTR model is inferred and overwrites the initial one, the option `marginal=True` can be used to construct a marginal rather than joint maximum likelihood reconstruction, and 'prune_short=False' can be used to avoid collapsing of zero length branches into polytomies.

  The tree and alignment arguments can be either file names (newick and fasta) or Biopython tree and alignent objects.

* Molecular clock phylogenies
    ```python
    from treetime import TreeTime
    tt = TreeTime(dates=mydates, 'my_tree.nwk', aln='my_seqs.nwk', gtr='JC69')
    tt.run(root='best', infer_gtr=True, relaxed_clock=(1.0,1.0), resolve_polytomies=True, max_iter=2)
    ```
  Every node of tt.tree will be assigned a `numdate` and `time_before_present` attribute. The additional attribute `resolve_polytomies` specifies whether TreeTime will attempt to resolve multiple mergers using the temporal constraints on leaves. Autocorrelated relaxed clocks can be fit by passing a tuple of two numbers `(slack, coupling)`. `slack` is the strength of the normal prior on rate variation, coupling penalizes rate variation between parents and children.

## Comparable Tools

There are several other tools which perform date inference for the internal nodes. [Beast](http://beast.bio.ed.ac.uk/) relies on the MCMC-type sampling of trees and is hence rather slow for large numbers of sequences. But BEAST allows the flexible inclusion of prior distributions, complex evolutionary models, and estimation of parameters.
Another tool that was recently published is [Least-Square-Dating](http://www.atgc-montpellier.fr/LSD/). LSD emphasises speed (also scales as O(N) as **TreeTime**), but provides limited scope for customization. **TreeTime** tries to achieve a good compromise betwee implementing a fast heuristic but accurate algorithm, which allows is readily extendable and customizable.


## Developer info

  - Credits -- .
  - Copyright and License: Pavel Sagulenko and Richard Neher, MIT Licence
  - How to contribute
  - References

