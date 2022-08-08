# TreeTime CLI overview

Many of TreeTime's commands define related workflows that contain each other partly.
The main difference between different workflows are whether they

* use sequence information (optional in some)
* use date information

The `homoplasy` and `mugration` are simple wrappers that can be added at a later point.

## `ancestral`

`ancestral` is a very fundamental workflow.
It requires sequence information and a tree and infers ancestral sequences for internal nodes on the tree and optionally
missing parts of sequence.

#### Inputs:

- ...

#### Steps:

- infer ancestral sequences
  - backward graph traversal, during which ...
  - calculate roots ...
  - forward graph traversal, during which ...

### Outputs:

- trees with internal nodes labeled
- inferred sequences

## `clock`

`clock` requires date information and in its basic usage does not require sequence information.
It only needs sequence information when run with the flag `--covariation`.
In this case, a more complicated algorithm is run that tries to estimate a clock model while accounting for shared
ancestry of samples.
In practice, the latter algorithm has not been very useful.

#### Inputs:

- ...

#### Steps:

- (optional): run clock-filter to label/exclude tips that are outliers
- (optional/only with `--covariation`): run a single iteration timetree analysis
- reroot the tree using `TreeTime.reroot`
- estimate clock model using `TreeTime.get_clock_model`

#### Outputs:

- rerooted tree
- table with dates and root-to-tip distances
- plot of root-to-tip distances vs time

## `timetree`

This is the default command run without subcommand.
It estimates a time tree and ultimately requires a tree (which can be inferred from an alignment).
It also requires either an alignment or the user to specify the sequence length.
This is needed for the probabilistic evolution models that determine how much a branch in the tree can be compressed and
stretched.
Much of the complexity of this workflow is hidden in the `TreeTime.run` function that goes through multiple iterations
of time tree estimation, rerooting, sequence inference, and so forth.
This complexity should be disentangled and made more explicit.

#### Inputs:

- ...

#### Steps

* `TreeTime.run` which under the hood does
  - initial iteration:
    - optimize the tree unless using the input branch length, otherwise infer ancestral sequences if an alignment is provided, remove unsupported branches.
    - (optional) reroot and clock filter. This is happening with the
    - `make_time_tree`: calculate branch length models either via inferring ancestral sequences (if sequences are provided) or using the branch length of the input tree and the sequence length
    - `make_time_tree`: infer initial time tree
  - subsequent iterations:
    - one initial reroot
    - optional step: coalescent models (most costly version only done in last iteration).
    - optional step: polytomy resolution
    - reinfer time tree and ancestral sequences
  - finally: evaluate rate variation and determine confidence intervals
* Happy path:
  - `infer_ancestral_sequences`
  - `reroot`
  - `make_time_tree`



#### Outputs:

* rerooted time tree
* plot of time tree
* table with dates and root-to-tip distances
* plot of root-to-tip distances vs time
* (optional): coalescent models
* molecular clock parameters

## `mugration`

This step is essentially a wrapper around `infer_ancestral_sequences` that pretends that traits like host species (pig,
chicken, etc) can be treated like sequence characters (A,C,G,T).
It assigns to each sample (tip in the tree) a sequence of length 1 with an alphabet that codes for the discrete traits.
These are often geographic locations and the name is a mash-up of "migration" and "mutation".

#### Inputs:

- ...

#### Steps

* coding of the discrete traits as one-character alphabet
* ancestral inference via the `reconstruct_discrete_traits` wrapper. This does some magic with GTR model inference that
  tries to counter sampling biases and basically inflates the rate at which characters change.

#### Outputs:

- ...

## `homoplasy`

Another wrapper around  `infer_ancestral_sequences`.
Here actual sequences are used and inferred.
The command then counts how often every site is mutated and compares this to a random distribution of mutations to
detect sites that are mutated suspiciously often.

#### Inputs:

- ...

#### Steps:

- ...

#### Outputs:

- ...

## `arg`

What is ARG?

#### Inputs:

- ...

#### Steps:

- ...

#### Outputs:

- ...
