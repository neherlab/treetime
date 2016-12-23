[![Build Status](https://travis-ci.org/neherlab/treetime.svg?branch=master)](https://travis-ci.org/neherlab/treetime)

# TreeTime: Maximum-Likelihood dating and ancestral inference for phylogenetic trees

## Overview

TreeTime provides routines for ancestral sequence reconstruction and the inference of molecular-clock phylogenies, i.e., a tree where all branches are scaled such that the locations of terminal nodes correspond to their sampling times and internal nodes are placed at the most likely time of divergence.

TreeTime aims at being a compromise between sophisticated probabilistic models of evolution and fast heuristics. It implements GTR models of ancestral inference and branch length optimization, but takes the tree topology as given.
The only topology optimization are resolution of polytomies in a way that is most (approximately) consistent with the sampling time constraints on the tree.
The package is designed to be used as a stand-alone tool or as a library used in larger phylogenetic analysis workflows.

#### Features
* ancestral sequence reconstruction (marginal and joint maximum likelihood)
* molecular clock tree inference (marginal and joint maximum likelihood)
* inference of GTR models
* rerooting to obtain best root-to-tip regression
* auto-correlated relaxed molecular clock (with normal prior)

![Molecular clock phylogeny of 200 NA sequences of influenza A H3N2](doc/flu_200.png)

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
TreeTime can be used as part of python programs that create and interact with tree time objects. The interface to treetime and its basic functionality are given below.

In addition, we provide scripts that can be run from the command line with arguments specifying input data and parameters.


* Ancestral sequence reconstruction:
  To perform ancestral sequence reconstruction, use the script `ancestral_inference.py`.
  ```
  usage: ancestral_reconstruction.py [-h] --aln ALN --tree TREE [--marginal]
                                     [--infer_gtr]

  Reconstructs ancestral sequences and maps mutations to the tree. The output
  consists of a file ending with _ancestral.fasta with ancestral sequences and a
  tree ending with _mutation.newick with mutations appended to node names as
  _A45G_... The inferred GTR model is written to stdout

  optional arguments:
    -h, --help   show this help message and exit
    --aln ALN    fasta file with input sequences
    --tree TREE  newick file with tree
    --marginal   marginal instead of joint ML reconstruction
    --infer_gtr  infer substitution model
  ```

  Alteratively, directly interact with the class TreeAnc from treetime as follows
    ```python
    from treetime import TreeAnc
    ta = TreeAnc(tree='my_tree.nwk', aln='my_seqs.nwk', gtr='JC69')
    ta.infer_ancestral_sequences(method = 'ml', infer_gtr=True, marginal=False)
    ```
  Every node of `ta.tree` now has a `node.sequence` attached. With the optional argument `infer_gtr=True`, a maximum likelihood GTR model is inferred and overwrites the initial one, the option `marginal=True` can be used to construct a marginal rather than joint maximum likelihood reconstruction, and 'prune_short=False' can be used to avoid collapsing of zero length branches into polytomies.

  The tree and alignment arguments can be either file names (newick and fasta) or Biopython tree and alignent objects.

* Molecular clock phylogenies
  To infer molecular clock phylogenies, use the script `timetree_inference.py`:
  ```
  usage: timetree_inference.py [-h] --aln ALN --tree TREE --dates DATES
                               [--infer_gtr] [--reroot REROOT]
                               [--resolve_polytomies]
                               [--relax [RELAX [RELAX ...]]]
                               [--max_iter MAX_ITER] [--verbose VERBOSE]
                               [--Tc TC] [--plot]

  Reconstructs ancestral sequences and infers a molecular clock tree. The output
  consists of a file ending with _ancestral.fasta with ancestral sequences and a
  tree ending with _mutation.newick with mutations appended to node names as
  _A45G_.... The branches of this tree are scaled such that branch length
  correspond to times in units of the molecular clock. The molecular clock,
  along with the inferred GTR model, is written to stdout

  optional arguments:
    -h, --help            show this help message and exit
    --aln ALN             fasta file with input sequences
    --tree TREE           newick file with tree
    --dates DATES         csv with dates for nodes with 'node_name, date' where
                          date is float (as in 2012.15)
    --infer_gtr           infer substitution model
    --reroot REROOT       reroot the tree. Valid arguments are 'best',
                          'midpoint', or a node name
    --resolve_polytomies  resolve polytomies using temporal information.
    --relax [RELAX [RELAX ...]]
                          use an autocorrelated molecular clock. Prior strength
                          and coupling of parent and offspring rates can be
                          specified e.g. as --relax 1.0 0.5
    --max_iter MAX_ITER   maximal number of iterations the inference cycle is
                          run
    --verbose VERBOSE     verbosity of output 0-6
    --Tc TC               coalescent time scale -- sensible values are on the
                          order of the average hamming distance of
                          contemporaneous sequences
    --plot                plot the tree on a time axis
  ```
  Alternatively, you can interact directly with the TreeTime class from within a script.

    ```python
    from treetime import TreeTime
    tt = TreeTime(dates=mydates, 'my_tree.nwk', aln='my_seqs.nwk', gtr='JC69')
    tt.run(root='best', infer_gtr=True, relaxed_clock=(1.0,1.0), resolve_polytomies=True, max_iter=2)
    ```
  Every node of tt.tree will be assigned a `numdate` and `time_before_present` attribute. The additional attribute `resolve_polytomies` specifies whether TreeTime will attempt to resolve multiple mergers using the temporal constraints on leaves. Autocorrelated relaxed clocks can be fit by passing a tuple of two numbers `(slack, coupling)`. `slack` is the strength of the normal prior on rate variation, coupling penalizes rate variation between parents and children.


* Quantify temporal signal in phylogenies and reroot to the maximize "clock-i-ness"
  The script `temporal_signal.py` provides functionality analogous to TempEst by Andrew Rambaut.
    ```
      usage: temporal_signal.py [-h] --tree TREE --dates DATES [--aln ALN]
                                [--infer_gtr] [--reroot] [--plot]
                                [--verbose VERBOSE]

      Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the
      tree. It will optionally reroot the tree to maximize the clock-like signal and
      recalculate branch length.

      optional arguments:
        -h, --help         show this help message and exit
        --tree TREE        newick file with tree
        --dates DATES      csv with dates for nodes with 'node_name, date' where
                           date is float (as in 2012.15)
        --aln ALN          fasta file with input sequences
        --infer_gtr        infer substitution model
        --reroot           reroot the tree to maximize the correlation of root-to-
                           tip distance with sampling time
        --plot
        --verbose VERBOSE  verbosity of output 0-6
    ```
  The slope of the regression of root-to-tip distance vs sampling date will be printed to stdout along with the fraction of variance explained by the linear regression. By passing the flag `--reroot`, treetime will search for the root that maximizes the correlation of root-to-tip distance with time and reroot the tree. The option `--plot` will produce a scatter plot with the best regression and save it to file.

## Comparable Tools

There are several other tools which estimate molecular clock phylogenies.
* [Beast](http://beast.bio.ed.ac.uk/) relies on the MCMC-type sampling of trees. It is hence rather slow for large data sets. But BEAST allows the flexible inclusion of prior distributions, complex evolutionary models, and estimation of parameters.
* [Least-Square-Dating](http://www.atgc-montpellier.fr/LSD/) (LSD) emphasises speed (it scales as O(N) as **TreeTime**), but provides limited scope for customization.


## Developer info

  - Credits -- .
  - Copyright and License: Pavel Sagulenko and Richard Neher, MIT Licence
  - How to contribute
  - References

