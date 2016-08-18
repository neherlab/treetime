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
        ta = TreeAnc(tree='my_tree.nwk', aln='my_seqs.nwk', gtr='Jukes-Cantor')
        ta.reconstruct_anc('ml')
    ```
Every node of `ta.tree` now has a `node.sequence` attached. Optimal arguments to 'reconstruct_anc' include `infer_gtr=True`, `marginal=True`, and 'prune_short=True'.



## Comparable Tools

There are several other tools which perform date inference for the internal nodes. [Beast](http://beast.bio.ed.ac.uk/) relies on the MCMC-type sampling of trees and is hence rather slow for large numbers of sequences. But BEAST allows the flexible inclusion of prior distributions, complex evolutionary models, and estimation of parameters.
Another tool that was recently published is [Least-Square-Dating](http://www.atgc-montpellier.fr/LSD/). LSD emphasises speed (also scales as O(N) as **TreeTime**), but provides limited scope for customization. **TreeTime** tries to achieve a good compromise betwee implementing a fast heuristic but accurate algorithm, which allows is readily extendable and customizable.


## Developer info

  - Credits -- .
  - Copyright and License: Pavel Sagulenko and Richard Neher, MIT Licence
  - How to contribute
  - References

