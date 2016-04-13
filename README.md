# TreeTime: Maximum-Likelihood dating and ancestral inference for phylogenetic trees

**treetime is currently in alpha status -- expect changes in the API**

## Overview

Estimate the date of the internal nodes of the phylogenetic trees using the maximum-likelihood algorithm given the temporal information about (some) leaves or internal nodes. The package takes the tree topology as an input in newick format, and a multiple-sequence alignment in fasta format. The temporal information about the leaves can be parsed from the fasta sequence annotations, or can be provided in a separate file. After parsing all the inputs, the dates for the rest of the nodes are inferred using maximum-likelihood algorithm in a manner similar to the ancestral state inference. Optionally, users can supply user-defined General-Time-Reversible models to perform custom maximum-likelihood optimization. Some of the standard models are supplied in the package and can give a good start of how to define a custom model.

In general, the algorithm relies on the given tree topology, and only optimizes the branch lengths with the exception of polytomies and zero-length branches. If the merging order is not resolved, the tree-building tools tend to create multiple mergers (so-called polytomies), or merge the unresolved branches in a random order. TreeTime reduces zero length branches to polytomies and resolves them in a way that is (approximately) most consistent with the sampling time constraints on the tree (trying to maximize the tree likelihood).

The package is designed to be used as a stand-alone tool as well as a module plugged in a bigger phylogenetic tool.


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

You might need root privileges for system wide installation. Alternatively, you can simply use it as-is. In this case, just download and unpack it, and then add the folder of the package location to your $PYTHONPATH variable.


### Basic usage


* Provide the function to extract the dates from the fasta annotation, or file with date-time information for the tree nodes. Say, we want to extract the dates directly from the alignment:
    ```python
    def fasta_name_to_date(fasta_name):
          # function definition
          pass
    ```

* Load the data from files to the TreeTime object:
```python
gtr = GTR.standard()
myTree = io.treetime_from_newick(gtr, nwk)
# set alignment to the tree
io.set_seqs_to_leaves(myTree, AlignIO.read(fasta, 'fasta'))
# set dates from the node names
io.set_node_dates_from_names(myTree, fasta_name_to_date) # note our custom date extractor here
```

* treetime requires a General time-reversible model (GTR) to estimate branch length and estimate location internal nodes. Very basic standard models are provided with gtr.py module, the above example uses a Jukes-Cantor model. If needed, the GTR class can be extended to supply a user-defined model of any complexity or estimated from the tree.

* Make unconstrained optimization  of the branch lengths and infer the ancestral sequences using ML algorithm (default). We perform iterative estimation of the ancestral sequences followed by optimization of the branch length. By default, This is especially important if there were (i) tree builder did not estimate the branch lengths consistently or (ii) a custom GTR model is used (distances might vary significantly).
    ```python
    myTree.optimize_seq_and_branch_len(infer_gtr=True)
    ```

* After the tree is fully prepared to the date inference, feed the date information to the leaves and run the optimization:
    ```python
    myTree.init_date_constraints()
    myTree.ml_t()
    ```

* Optionally, resolve polytomies:
    ```python
    myTree.resolve_polytomies()
    ```

* Save the tree:
    ```python
    io.save_tree(myTree, format='newick', write_dates=True)
    ```


The full example for the basic usage of the package can be found in the examples section.


## Design Goals

The purpose for this package is to provide better, more realistic phylogenetic trees. It is intended to work in pair with the standard tools for the phylogenetic unferrence (as a sub-module, or stand-alone). We targeted the robust and fast tool for general purpose. Given these goals, we aim to create the tool which is

(i) fast and lightweight. Currently, the algorithm scales with the number of tree leaves as O(N).

(ii) flexible enough to be plugged into different bigger software modules and programms.

(iii) We want to make a tool, which can be easily extended to any level of complexity basing on particular requirements.


## Comparable Tools

There are several other tools which perform date inference for the internal nodes. [Beast](http://beast.bio.ed.ac.uk/) relies on the MCMC-type sampling of trees and is hence rather slow for large numbers of sequences. But BEAST allows the flexible inclusion of prior distributions, complex evolutionary models, and estimation of parameters.
Another tool that was recently published is [Least-Square-Dating](http://www.atgc-montpellier.fr/LSD/). LSD emphasises speed (also scales as O(N) as **TreeTime**), but provides limited scope for customization. **TreeTime** tries to achieve a good compromise betwee implementing a fast heuristic but accurate algorithm, which allows is readily extendable and customizable.


## Developer info

  - Credits -- .
  - Copyright and License: Pavel Sagulenko and Richard Neher, MIT Licence
  - How to contribute
  - References

