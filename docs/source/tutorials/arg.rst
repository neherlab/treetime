
Using Recombination Event Knowledge to Improve Time Tree Inference
------------------------------------------------------------------

Although relatively uncommon, recombination can be a driver of pathogen evolution. However, recombining segments are often not used together in phylogenetic inference due to computational challenges.
Because of this, typically, segments of the genome with little to no recombination are used to infer the pathogen phylogeny, leading to a loss of information. 
As recombination is an uncommon event most segments will have portions of their phylogenies that show a great deal of overlap. Segments share evolutionary history in these areas of overlap and this knowledge can 
be used to improve divergence time estimates in TreeTime.

To infer branch length TreeTime makes the assumption that the number of mutations on a branch of length :math:`\tau` is Poisson distributed 

.. math::
    n_{mut} \sim Pois(\mu \tau). 

Where :math:`\mu` is the mutation rate. 
The variance and mean number of expected mutations on a branch of length :math:`\tau` is :math:`\mu \tau`.
If :math:`\mu` is known this relation can be used to estimate the branch length given the number of seen mutations, 

.. math::
    \tau^{infer} = \frac{n_{mut}}{\mu},

this estimator has expectation :math:`\tau`` and variance :math:`\frac{\tau}{\mu}`.
When this is repeated :math:`L` times (i.e. for a sequence of :math:`L` nucleotides, each with mutation rate :math:`\mu`), 
the average number of mutations is normally distributed with mean :math:`\mu \tau` and standard error :math:`\frac{\mu \tau}{\sqrt{L}}`. 
Assuming we know that a branch is shared between two segments we can use the alignment of both segments on this branch to estimate divergence time of this branch, decreasing the standard error.

`TreeKnit <https://github.com/PierreBarrat/TreeKnit.jl>`_ is a package that can infer recombination events from tree topologies. It returns lists of leaves that are connected by shared branches in pairs of trees. 
These leaves can be used to determine so called maximally compatible clades (MCCs), or clades where topology is shared across trees.  
If desired TreeKnit additionally returns trees that have been resolved according to each other. 

This output can be used in TreeTime to improve the inference of time trees and in turn improve the ancestral sequence reconstruction and the clock tree inference. 
The ``treetime arg`` command uses input trees and their corresponding alignments to infer time trees. It is assumed the list of MCCs are in json format as described in TreeKnit. 
For each tree, the list of maximally compatible clades with every other tree is used to determine if internal nodes are part of a MCC with another tree and if they are, which MCC they belong to. 
This is done using the Fitch algorithm (function: ``assign_all_mccs``). If a node and it's parent both belong to the same MCC then the branch between them is shared.
For example if a branch is shared between trees of segments A and B then the alignment of both segment A and B can be used to infer the divergence time of the branch, 
leading to more accurate branch length estimates than if only the alignment of segment A was used to infer the divergence time of that branch.

In the test folder there is an example of a standard TreeKnit output for three trees. TreeTime expects recombination information to be in `TreeKnit output format <https://pierrebarrat.github.io/TreeKnit.jl/overview/#Output>`_. This can be used to run the ``treetime arg`` command:

.. code:: bash

    treetime arg --trees arg/TreeKnit/tree_a_resolved.nwk arg/TreeKnit/tree_b_resolved.nwk arg/TreeKnit/tree_c_resolved.nwk --alignments arg/TreeKnit/aln_a.fasta arg/TreeKnit/aln_b.fasta arg/TreeKnit/aln_c.fasta --mccs arg/TreeKnit/MCCs.json --dates arg/TreeKnit/metadata.csv --clock-rate 0.0028 --outdir time_tree_arg_results

For each tree treetime will output ancestral sequence reconstructions, dates of the tree nodes, as well as time tree and divergence trees for each input tree using information from other trees for shared branches. 
The output will be written to the folder ``time_tree_arg_results``.
