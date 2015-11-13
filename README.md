# TreeTime: Maximum-Likelihood datign of the internal nodes in phylogenetic trees 


## Overview

Estimate the date of the internal nodes of the phylogenetic trees using the maximum-likelihood algorithm given the temporal information about (some) leaves or internal nodes. The package takes the tree topology as an input in newick format, and a multiple-sequence alignment in fasta format. The temporal information about the leaves can be parsed from the fasta sequence annotations, or can be provided in a separate file. After parsing all the inputs, the dates for the rest of the nodes are inferred using maximum-likelihood algorithm in a manner similar to the ancestral state inferrence. Optionally, users can supply user-defined General-Time-Reversible models to perform custom maximum-likelihood optimization. Some of the standard models are supplied in the package and can give a good start of how to define a custom model. 

In general, the algorithm relies on teh given tree topology, and only optimizes the branch lengths. There is however, one exception. If the merging order is not resolved, the tree-building tools tend to create multiple mergers (so-called polytomies), or merge the unresolved branches in a reandom order. Our tool is able to resolve the polytomies to merge the nodes in the order which maximizes the tree likelihood. 

The package is desinged to be used as a standalone tool as well as a module plugged in a bigger phylogenetic tool. 


## Getting started

### Instalation and prerequisites

* The package depends on several python libraries: 
    - numpy, SciPy: for all kind of mathematical operations as matrix operations, numerical integration, interpolation, minimization, etc. 
  
    - BioPython: for multiple sequence alignment and all phylogenetic functionality
  If you do not have these libraries, you can install them by typing in the terminal: 
    ```bash
    $pip install numpy scipy biopython
    ```

* To install the package, run install.sh script from the terminal: 
    ```bash
    $./install.sh
    ```

Or, you can simply use it as-is. In this case, just download and unpack it, and then add the folder of the package location to your $PYTHONPATH variable.


### Basic usage


* Provide the function to extract the dates from the fasta annotation, or file with date-time information for the tree nodes. Say, we want to extract the dates directly from the alinment:
    ```python
    def fasta_name_to_date(fasta_name):
          # function definition
          pass
    ```

* Load the data from files to the TreeTime object:
```python
tree_time = io.treetime_from_newick(nwk)
# set alignment to the tree
io.set_seqs_to_leaves(tree_time, AlignIO.read(fasta, 'fasta'))
# set dates from the node names
io.set_node_dates_from_names(tree_time, fasta_name_to_date) # note our custom date extractor here
```

* Then create the General time-reversible model (GTR) object. There are some provided with gtr.py module, we gonna use the default one (Jukes-Cantor model). If needed, the GTR class can be extended to supply a user-defined model of any complexity. 
    ```python
    gtr = GTR.standard()    
    ```

* Make unconstrained optimization  of the branch lengths and infer the ancestral sequences using ML algortihm (default). If you rely on the tree quality (mostly, on the branch lengths), the tree branch lenghts optimization can be skipped. Despite it provides some performance benefits, it is however discouraged to do so because this step binds the units of time used by the GR model and the tree raw branch lengths. This is especially important if there were (i) tree builder did not estimate the branch lengths consistently or (ii) a custom GTR model is used (distances might vary significantly). 
    ```python
    tree_time.optimize_seq_and_branch_len(gtr)
    ```

* After the tree is fully prepared to the date inferrence, feed the date information to the leaves and run the optimization:
    ```python
    tree_time.init_date_constraints(gtr, slope=slope)
    tree_time.ml_t(gtr)
    ```

* Optionally, resolve polytomies: 
    ```python
    tree_time.resolve_polytomies(gtr)
    ```

* Save the tree:
    ```python
    io.save_tree(tree_time, format='newick', write_dates=True)
    ```


The full example for the basic usage of the package can be found in the examples section. 


## Design Goals

The purpose for this package is to provide better, more realistic phylogenetic trees. It is intended to work in pair with the standard tools for the phylogenetic unferrence (as a sub-module, or stand-alone). We targeted the robust and fast tool for general purpose. Given these goals, we aim to create the tool which is 

(i) fast and lightweight. Currently, the algorithm scales with the number of tree leaves as O(N). 

(ii) flexible enough to be plugged into different bigger software modules and programms. 

(iii) We want to make a tool, which can be easily extended to any level of complexity basing on particular requirements. 


## Comparable Tools

There are few other tools which perform date inferrence for the internal nodes. Worth noting [Beast](http://beast.bio.ed.ac.uk/). This tool relies on the MCMC-type algorithm, which makes its rather slow, though it allows for bigger functionality and more thorough results. 
Another tool is recently appeared [Least-Square-Dating](http://www.atgc-montpellier.fr/LSD/), which is the opposite to BEAST: it is extremely fast (also scales as O(N) as **TreeTime**), but is poorly extendable. We provide the trade-off between the two implementing the fast algorithm, which allows for any type of extending.


## Developer info

  - Credits -- .
  - Copyright and License -- 
  - How to contribute
  - References

