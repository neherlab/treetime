---
sidebar_position: 9999
---

# Command-line reference

<!-- WARNING: GENERATED CONTENT -->
<!-- Do not edit the generated file manually! All manual changes will be overwritten by automation. -->
<!-- See developer guide for how to properly work with this. -->

:::warning
Generated content. Report issues at https://github.com/neherlab/treetime/issues
:::




This document contains the automatically generated reference documentation for command-line arguments of the latest version of TreeTime CLI.

If you have TreeTime CLI installed, you can type `treetime --help` to read the latest documentation for your installed version of TreeTime. To generate this document in markdown format, run `treetime help-markdown > reference.md`
  

**Command Overview:**

* [`treetime`↴](#treetime)
* [`treetime completions`↴](#treetime-completions)
* [`treetime timetree`↴](#treetime-timetree)
* [`treetime optimize`↴](#treetime-optimize)
* [`treetime prune`↴](#treetime-prune)
* [`treetime ancestral`↴](#treetime-ancestral)
* [`treetime clock`↴](#treetime-clock)
* [`treetime homoplasy`↴](#treetime-homoplasy)
* [`treetime mugration`↴](#treetime-mugration)
* [`treetime arg`↴](#treetime-arg)
* [`treetime schema`↴](#treetime-schema)
* [`treetime help-markdown`↴](#treetime-help-markdown)

## `treetime`

Maximum-likelihood phylodynamic inference

Documentation: https://treetime.readthedocs.io/en/stable/
Publication:   https://academic.oup.com/ve/article/4/1/vex042/4794731

**Usage:** `treetime [OPTIONS] <COMMAND>`

###### **Subcommands:**

* `completions` — Generate shell completions
* `timetree` — Estimates time trees from an initial tree topology, a set of date constraints (e.g. tip dates), and an alignment (optional)
* `optimize` — Optimizes the branch lengths and likelihood of a phylogenetic tree given aligned sequences
* `prune` — Prunes short branches and/or branches without mutations from a phylogenetic tree
* `ancestral` — Reconstructs ancestral sequences and maps mutations to the tree. The output consists of a file 'ancestral.fasta' with ancestral sequences and a tree 'annotated_tree.nexus' with mutations added as comments like A45G,G136T,..., number in SNPs used 1-based index by default. The inferred GTR model is written to stdout
* `clock` — Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. It will reroot the tree to maximize the clock-like signal and recalculate branch length unless run with --keep_root
* `homoplasy` — Reconstructs ancestral sequences and maps mutations to the tree. The tree is then scanned for homoplasies. An excess number of homoplasies might suggest contamination, recombination, culture adaptation or similar
* `mugration` — Reconstructs discrete ancestral states, for example geographic location, host, or similar. In addition to ancestral states, a GTR model of state transitions is inferred
* `arg` — Estimates ancestral reassortment graph (ARG)
* `schema` — Write JSON Schema definitions for TreeTime data types
* `help-markdown` — Print CLI reference documentation in Markdown format

###### **Options:**

* `-j`, `--jobs <JOBS>` — Number of processing jobs. If not specified, all available CPU threads will be used

  Default value: `20`
* `--verbosity <VERBOSITY>` — Set verbosity level of console output

  Default value: `warn`

  Possible values: `off`, `error`, `warn`, `info`, `debug`, `trace`

* `--silent` — Disable all console output. Same as `--verbosity=off`
* `-v`, `--verbose` — Make console output more verbose. Add multiple occurrences to increase verbosity further
* `-q`, `--quiet` — Make console output more quiet. Add multiple occurrences to make output even more quiet
* `--no-progress` — Disable progress bar display



## `treetime completions`

Generate shell completions.

This will print the completions file contents to the console. Refer to your shell's documentation on how to install the completions.

Example for Ubuntu Linux:

treetime completions bash > ~/.local/share/bash-completion/treetime

**Usage:** `treetime completions [SHELL]`

###### **Arguments:**

* `<SHELL>` — Name of the shell to generate appropriate completions

  Default value: `bash`

  Possible values: `bash`, `elvish`, `fish`, `fig`, `powershell`, `zsh`




## `treetime timetree`

Estimates time trees from an initial tree topology, a set of date constraints (e.g. tip dates), and an alignment (optional)

**Usage:** `treetime timetree [OPTIONS] --outdir <OUTDIR>`

###### **Options:**

* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Path to one or multiple FASTA files with aligned input sequences

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

   See: https://en.wikipedia.org/wiki/FASTA_format
* `-t`, `--tree <TREE>` — Name of file containing the tree in newick, nexus, or phylip format.

   If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
* `-r`, `--vcf-reference <VCF_REFERENCE>` — Only for vcf input: fasta file of the sequence the VCF was mapped to
* `-d`, `--metadata <METADATA>` [alias: `dates`] — CSV/TSV file with metadata including sampling dates
* `--metadata-id-columns <COLUMN>` [alias: `name-column`] — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma and tab.

  Default values: `,`, `	`
* `--date-column <COLUMN>` — Label of the column to be used as sampling date (auto-detected when omitted)
* `--date-format <FORMAT>` — Format used to parse string sampling dates (numeric, ISO, and uncertain dates parse regardless)

  Default value: `%Y-%m-%d`
* `--sequence-length <SEQUENCE_LENGTH>` — Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided
* `--clock-rate <CLOCK_RATE>` — If specified, the rate of the molecular clock won't be optimized
* `--clock-std-dev <CLOCK_STD_DEV>` — Standard deviation of the provided clock rate estimate
* `--branch-length-mode <BRANCH_LENGTH_MODE>` — If set to 'input', the provided branch length will be used without modification. Branch lengths optimized by treetime are only accurate at short evolutionary distances

  Default value: `marginal`

  Possible values: `input`, `marginal`

* `--time-marginal <TIME_MARGINAL>` — Control when marginal time distributions are used for output.

   All modes use marginal inference during optimization. The mode controls whether confidence intervals are extracted from the resulting distributions:

   - `never`: no confidence interval output (default) - `always`: write confidence intervals from distributions computed during optimization - `only-final`: run one extra inference pass after optimization, then write confidence intervals

  Default value: `never`

  Possible values: `never`, `always`, `only-final`

* `--confidence` — Add rate-uncertainty to confidence intervals.

   `--time-marginal=always` and `only-final` already write mutation-stochasticity CIs. This flag adds rate-uncertainty CIs (re-runs inference at rate +/- sigma), combined via quadrature sum. Requires `--covariation` or `--clock-std-dev`. When set with `--time-marginal=never` (default), automatically promotes to `only-final`.
* `--keep-polytomies` — Don't resolve polytomies using temporal information
* `--resolve-polytomies` — Resolve polytomies using temporal information
* `--relax <SLACK>` — use an autocorrelated molecular clock. Strength of the gaussian priors on branch specific rate deviation and the coupling of parent and offspring rates can be specified e.g. as --relax 1.0 0.5. Values around 1.0 correspond to weak priors, larger values constrain rate deviations more strongly. Coupling 0 (--relax 1.0 0) corresponds to an un-correlated clock
* `--max-iter <MAX_ITER>` — maximal number of iterations the inference cycle is run. For polytomy resolution and coalescence models max_iter should be at least 2

  Default value: `2`
* `--coalescent <COALESCENT>` — Coalescent time scale in years.

   Sensible values are on the order of the average hamming distance of contemporaneous sequences divided by the clock rate. For example, if average pairwise distance is 0.01 substitutions/site and clock rate is 0.001 subs/site/year, then Tc ~ 10 years.
* `--coalescent-opt` — Optimize coalescent time scale Tc to maximize coalescent likelihood.

   When set, TreeTime will find the optimal Tc using Brent's method. This is equivalent to Python v0's `--coalescent=opt`. If --coalescent is also provided, that value is used as the initial guess; otherwise defaults to 1.0.
* `--coalescent-skyline` — Use skyline coalescent model instead of constant Tc.

   Estimates a piecewise linear coalescent rate history. Requires --n-skyline to specify the number of grid points.
* `--n-skyline <N_SKYLINE>` — Number of grid points in skyline coalescent model.

   Only used when --coalescent-skyline is set. Defines how many piecewise linear segments are used to model Tc(t) over time. Must be at least 2.

  Default value: `10`
* `--tip-labels` — add tip labels (default for small trees with <30 leaves)
* `--no-tip-labels` — don't show tip labels (default for trees with >=30 leaves)
* `--clock-filter <CLOCK_FILTER>` — ignore tips that don't follow a loose clock, 'clock-filter=number of inter-quartile ranges from regression'. Default=3.0, set to 0 to switch off

  Default value: `3.0`
* `--n-iqd <N_IQD>` — Number of IQD (interquartile distance) for clock filter outlier detection
* `--reroot <REROOT>` — Reroot the tree using root-to-tip regression. Valid choices are 'min_dev', 'least-squares', and 'oldest'. 'least-squares' adjusts the root to minimize residuals of the root-to-tip vs sampling time regression, 'min_dev' minimizes variance of root-to-tip distances. 'least- squares' can be combined with --covariation to account for shared ancestry. Alternatively, you can specify a node name or a list of node names to be used as outgroup or use 'oldest' to reroot to the oldest node. By default, TreeTime will reroot using 'least-squares'. Use --keep- root to keep the current root

  Default value: `least-squares`

  Possible values: `least-squares`, `min-dev`, `oldest`, `clock-filter`, `mrca`

* `--keep-root` — don't reroot the tree. Otherwise, reroot to minimize the residual of the regression of root-to-tip distance and sampling time
* `--allow-negative-rate` — By default, rates are forced to be positive. For trees with little temporal signal it is advisable to remove this restriction to achieve essentially mid-point rooting
* `--tip-slack <TIP_SLACK>` — excess variance associated with terminal nodes accounting for overdispersion of the molecular clock
* `--covariation` — Account for covariation when estimating rates or rerooting using root-to-tip regression
* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `--model infer` infers a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use `--model-params` to specify those.

  Default value: `infer`

  Possible values:
  - `infer`:
    Infer GTR parameters from data via Fitch parsimony substitution counts
  - `jc69`
  - `k80`
  - `f81`
  - `hky85`
  - `t92`
  - `tn93`
  - `jtt92`

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as a `key=value` list

   Example: `--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.

   See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`, `joint`

* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--dense <DENSE>` — Use dense representation for sequences (store full probability distributions)

  Possible values: `true`, `false`

* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--zero-based` — Zero-based mutation indexing
* `--reconstruct-tip-states` — Overwrite ambiguous states on tips with the most likely inferred state
* `--report-ambiguous` — Include transitions involving ambiguous states
* `--no-indels` — Disable indel (insertion/deletion) contributions to branch-length optimization and branch-length distributions.

   When set, branch-length optimization uses substitution-only likelihood and timetree branch distributions exclude the Poisson indel term. Matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and enables v0 parity testing. Default: indels enabled.
* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Write augur-compatible node data JSON to this path

   Contains per-node dates, branch lengths, clock model parameters, confidence intervals, and divergence metrics. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration. Defaults to `<outdir>/timetree.augur-node-data.json`.
* `-O`, `--outdir <OUTDIR>` — Directory to write the standard set of output files to
* `--tracelog <TRACELOG>` — Write iteration statistics to tracelog CSV file for monitoring convergence
* `--seed <SEED>` — Random seed



## `treetime optimize`

Optimizes the branch lengths and likelihood of a phylogenetic tree given aligned sequences

**Usage:** `treetime optimize [OPTIONS] --tree <TREE> --outdir <OUTDIR>`

###### **Options:**

* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Path to one or multiple FASTA files with aligned input sequences

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

   See: https://en.wikipedia.org/wiki/FASTA_format
* `-t`, `--tree <TREE>` — Name of file containing the tree in newick, nexus, or phylip format.

   If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `--model infer` infers a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use `--model-params` to specify those.

  Default value: `infer`

  Possible values:
  - `infer`:
    Infer GTR parameters from data via Fitch parsimony substitution counts
  - `jc69`
  - `k80`
  - `f81`
  - `hky85`
  - `t92`
  - `tn93`
  - `jtt92`

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as a `key=value` list

   Example: `--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.

   See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
* `--dense <DENSE>` — Use dense representation of sequences on the tree

   Dense mode stores full probability vectors at every alignment position for each node. Sparse mode stores only variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

  Possible values: `true`, `false`

* `-O`, `--outdir <OUTDIR>` — Directory to write the standard set of output files to
* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Write augur-compatible node data JSON to this path

   Contains per-node optimized branch lengths (divergence, substitutions per site) and the input alignment and tree paths. The output is compatible with augur export v2 --node-data, equivalent to `augur refine` run without `--timetree` (no clock or date fields). Defaults to `<outdir>/optimize.augur-node-data.json`.
* `--max-iter <MAX_ITER>` — Maximum number of iterations

  Default value: `10`
* `--dp <DP>` — Likelihood convergence threshold. The loop stops when successive likelihoods differ by less than this value, or when a 2-cycle with amplitude below this value is detected

  Default value: `0.1`
* `--damping <DAMPING>` — Damping factor for outer-loop branch length updates.

   Controls how aggressively new branch lengths replace old ones during iterative optimization. At each iteration i, the update is: bl = bl_new * (1 - d) + bl_old * d where d = max(damping^(i+1), 0.01). The 1% floor prevents fully undamped late iterations on non-monotone objectives.

   Higher values are more conservative (slower convergence, less oscillation). Set to 0.0 to disable damping (full update each iteration, no floor). Must be in [0.0, 1.0).

  Default value: `0.75`
* `--branch-length-initial-guess <BRANCH_LENGTH_INITIAL_GUESS>` — Initial branch length estimate before Newton optimization.

   - auto: estimate only edges with missing or invalid branch lengths, preserve valid input values (default) - always: estimate all edges, overwriting input branch lengths - never: use input branch lengths as-is; fails if any are missing

  Default value: `auto`

  Possible values:
  - `auto`:
    Estimate only edges with missing or invalid branch lengths, preserve valid input values. No-op when all edges have finite branch lengths
  - `always`:
    Estimate all edges, overwriting input branch lengths
  - `never`:
    Use input branch lengths as-is. Fails if any edge has a missing or invalid branch length

* `--opt-method <OPT_METHOD>` — Per-edge branch length optimization method.

   Algorithm x parameterization: - brent: Brent's method in t space (derivative-free) - brent-sqrt: Brent's method in sqrt(t) space (default, matches v0) - brent-log: Brent's method in ln(t) space - newton: Newton-Raphson in t space - newton-sqrt: Newton-Raphson in sqrt(t) space - newton-log: Newton-Raphson in ln(t) space

  Default value: `brent-sqrt`

  Possible values:
  - `brent`:
    Brent's method in $t$ space (derivative-free, bracket-based)
  - `brent-sqrt`:
    Brent's method in $\sqrt{t}$ space
  - `brent-log`:
    Brent's method in $\ln(t)$ space
  - `newton`:
    Newton-Raphson in $t$ space
  - `newton-sqrt`:
    Newton-Raphson in $\sqrt{t}$ space
  - `newton-log`:
    Newton-Raphson in $\ln(t)$ space

* `--no-indels` — Disable indel (insertion/deletion) contributions to branch-length optimization.

   When set, the optimizer uses substitution-only likelihood, matching standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and enabling v0 parity testing. Default: indels enabled.
* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`




## `treetime prune`

Prunes short branches and/or branches without mutations from a phylogenetic tree

**Usage:** `treetime prune [OPTIONS] --tree <TREE> --outdir <OUTDIR>`

###### **Options:**

* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Path to one or multiple FASTA files with aligned input sequences

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

   See: https://en.wikipedia.org/wiki/FASTA_format
* `-t`, `--tree <TREE>` — Name of file containing the tree in newick, nexus, or phylip format
* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `-O`, `--outdir <OUTDIR>` — Directory to write the standard set of output files to
* `-s`, `--prune-short <THRESHOLD>` — Threshold value for pruning of branches

   If set, prune branches with a length below this value
* `-e`, `--prune-empty` — Prune empty branches

   If set, prune any branch that does not have a mutation or other state transition mapped to it.

   Requires --alignment
* `-m`, `--merge-shared-mutations` — Merge branches in polytomies that share identical mutations

   When sibling branches in a polytomy (node with >2 children) carry identical substitutions, they are grouped under a new internal node. The shared mutations move to the new branch (parent to new node), and only unique mutations remain on children's edges. Reduces tree builder artifacts from arbitrary binary resolution of polytomies.

   Requires --alignment
* `-n`, `--prune-nodes-list <NODE_NAMES>` — List of node names to prune

   List of node names to remove from the tree, comma-separated (,)

   Use --prune-nodes-list-delimiter to specify a different delimiter.
* `--prune-nodes-list-delimiter <DELIMITER>` — Name separator for `--prune-nodes-list`

   String used to separate node names in the list given to (--prune-nodes-list). Make sure to correctly quote and escape the delimiter according to your shell.

  Default value: `,`
* `-N`, `--prune-nodes-list-file <FILEPATH>` — File containing list of node names to prune

   Path to a file containing node names to remove from the tree, newline-delimited (\n).

   Use '-' to read from standard input (stdin).

   Use --prune-nodes-list-file-delimiter to specify a different delimiter.
* `--prune-nodes-list-file-delimiter <DELIMITER>` — Separator for node names in the list file

   Character or string used to separate node names in the list file.

  Default value: `
`



## `treetime ancestral`

Reconstructs ancestral sequences and maps mutations to the tree. The output consists of a file 'ancestral.fasta' with ancestral sequences and a tree 'annotated_tree.nexus' with mutations added as comments like A45G,G136T,..., number in SNPs used 1-based index by default. The inferred GTR model is written to stdout

**Usage:** `treetime ancestral [OPTIONS] --tree <TREE> --outdir <OUTDIR>`

###### **Options:**

* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Path to one or multiple FASTA files with aligned input sequences

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

   See: https://en.wikipedia.org/wiki/FASTA_format
* `-r`, `--vcf-reference <VCF_REFERENCE>` — FASTA file of the sequence the VCF was mapped to (only for vcf input)
* `-t`, `--tree <TREE>` — Name of file containing the tree in newick, nexus, or phylip format.

   If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `--model infer` infers a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use `--model-params` to specify those.

  Default value: `infer`

  Possible values:
  - `infer`:
    Infer GTR parameters from data via Fitch parsimony substitution counts
  - `jc69`
  - `k80`
  - `f81`
  - `hky85`
  - `t92`
  - `tn93`
  - `jtt92`

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as a `key=value` list

   Example: `--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.

   See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`, `joint`

* `--dense <DENSE>` — Use dense representation (stores full probability vectors at each position)

   When combined with `--model infer`, marginal reconstruction runs twice: once to populate profiles for GTR inference, and again with the inferred GTR.

  Possible values: `true`, `false`

* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--zero-based` — Zero-based mutation indexing
* `--reconstruct-tip-states` — Overwrite ambiguous states on tips with the most likely inferred state
* `--report-ambiguous` — Include transitions involving ambiguous states
* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Write augur-compatible node data JSON to this path.

   Contains per-node nucleotide mutations, reconstructed sequences, the alignment mask, genome annotations, and the reference (root) sequence. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.
* `-O`, `--outdir <OUTDIR>` — Directory to write the standard set of output files to
* `--gtr-iterations <GTR_ITERATIONS>` — Number of outer GTR refinement iterations.

   Re-estimates the rate matrix from marginal posterior profiles after each reconstruction pass. Only effective with `--model infer`. Default 0 preserves the current single-pass behavior. Mugration uses 5 by default.

  Default value: `0`
* `--seed <SEED>` — Random seed
* `--sample-from-profile <SAMPLE_FROM_PROFILE>` — How to pick ancestral states from the marginal posterior profile.

   'argmax': most likely state at every node (deterministic, default). 'root': sample from the posterior at the root only, argmax elsewhere (matches augur's `sample_from_profile='root'`). Use `--seed` for reproducible draws. 'all': sample from the posterior at every node.

   Only affects marginal reconstruction (`--method-anc=marginal`).

  Default value: `argmax`

  Possible values: `argmax`, `root`, `all`




## `treetime clock`

Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. It will reroot the tree to maximize the clock-like signal and recalculate branch length unless run with --keep_root

**Usage:** `treetime clock [OPTIONS] --metadata <METADATA> --outdir <OUTDIR>`

###### **Options:**

* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Path to one or multiple FASTA files with aligned input sequences

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

   See: https://en.wikipedia.org/wiki/FASTA_format
* `-t`, `--tree <TREE>` — Name of file containing the tree in newick, nexus, or phylip format.

   If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
* `-r`, `--vcf-reference <VCF_REFERENCE>` — Only for vcf input: fasta file of the sequence the VCF was mapped to
* `-d`, `--metadata <METADATA>` [alias: `dates`] — CSV/TSV file with metadata including sampling dates
* `--metadata-id-columns <COLUMN>` [alias: `name-column`] — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma and tab.

  Default values: `,`, `	`
* `--date-column <COLUMN>` — Label of the column to be used as sampling date (auto-detected when omitted)
* `--date-format <FORMAT>` — Format used to parse string sampling dates (numeric, ISO, and uncertain dates parse regardless)

  Default value: `%Y-%m-%d`
* `--sequence-length <SEQUENCE_LENGTH>` — Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided
* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `--model infer` infers a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use `--model-params` to specify those.

  Default value: `infer`

  Possible values:
  - `infer`:
    Infer GTR parameters from data via Fitch parsimony substitution counts
  - `jc69`
  - `k80`
  - `f81`
  - `hky85`
  - `t92`
  - `tn93`
  - `jtt92`

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as a `key=value` list

   Example: `--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.

   See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
* `--branch-length-mode <BRANCH_LENGTH_MODE>` — If set to 'input', the provided branch length will be used without modification. Note that branch lengths optimized by treetime are only accurate at short evolutionary distances

  Default value: `marginal`

  Possible values: `input`, `marginal`

* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`, `joint`

* `--clock-filter <CLOCK_FILTER>` — ignore tips that don't follow a loose clock, 'clock-filter=number of interquartile ranges from regression'. Default=3.0, set to 0 to switch off

  Default value: `3.0`
* `--reroot <REROOT>` — Reroot the tree using root-to-tip regression. Valid choices are 'min_dev', 'least-squares', and 'oldest'. 'least-squares' adjusts the root to minimize residuals of the root-to-tip vs sampling time regression, 'min_dev' minimizes variance of root-to-tip distances. 'least- squares' can be combined with --covariation to account for shared ancestry. Alternatively, you can specify a node name or a list of node names to be used as outgroup or use 'oldest' to reroot to the oldest node. By default, TreeTime will reroot using 'least-squares'. Use --keep- root to keep the current root

  Default value: `least-squares`

  Possible values: `least-squares`, `min-dev`, `oldest`, `clock-filter`, `mrca`

* `--keep-root` — don't reroot the tree. Otherwise, reroot to minimize the residual of the regression of root-to-tip distance and sampling time
* `--prune-short`
* `--tip-slack <TIP_SLACK>` — excess variance associated with terminal nodes accounting for overdispersion of the molecular clock
* `--covariation` — Account for covariation when estimating rates or rerooting using root-to-tip regression
* `--allow-negative-rate` — By default, rates are forced to be positive. For trees with little temporal signal it is advisable to remove this restriction to achieve essentially mid-point rooting
* `-O`, `--outdir <OUTDIR>` — Directory to write the standard set of output files to
* `--seed <SEED>` — Random seed
* `--branch-split-method <METHOD>` — Optimization method to use for finding the best root position

  Default value: `grid`

  Possible values:
  - `grid`:
    Grid search with equally-spaced evaluation points
  - `brent`:
    Brent's method for robust 1D optimization
  - `golden-section`:
    Golden section search optimization

* `--branch-split-grid-n-points <N_POINTS>` — Number of equally-spaced points to evaluate (grid method only)

  Default value: `11`
* `--branch-split-brent-max-iters <BRENT_MAX_ITERS>` — Maximum number of iterations for Brent's method

  Default value: `50`
* `--branch-split-brent-tolerance <BRENT_TOLERANCE>` — Convergence tolerance for Brent's method

  Default value: `0.000000000001`
* `--branch-split-golden-max-iters <GOLDEN_MAX_ITERS>` — Maximum number of iterations for golden section search

  Default value: `50`
* `--branch-split-golden-tolerance <GOLDEN_TOLERANCE>` — Convergence tolerance for golden section search

  Default value: `0.000000000001`
* `--variance-factor <VARIANCE_FACTOR>` — Variance scaling factor proportional to branch length

  Default value: `0`
* `--variance-offset <VARIANCE_OFFSET>` — Constant variance offset for all branches

  Default value: `0`
* `--variance-offset-leaf <VARIANCE_OFFSET_LEAF>` — Additional variance offset for leaf (terminal) nodes

  Default value: `1`



## `treetime homoplasy`

Reconstructs ancestral sequences and maps mutations to the tree. The tree is then scanned for homoplasies. An excess number of homoplasies might suggest contamination, recombination, culture adaptation or similar

**Usage:** `treetime homoplasy [OPTIONS] --tree <TREE> --outdir <OUTDIR>`

###### **Options:**

* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Path to one or multiple FASTA files with aligned input sequences

   Accepts plain or compressed FASTA files. If a compressed fasta file is provided, it will be transparently decompressed. Supported compression formats: `gz`, `bz2`, `xz`, `zstd`. Decompressor is chosen based on file extension. If there's multiple input files, then different files can have different compression formats.

   If no input files provided, the plain fasta input is read from standard input (stdin).

   See: https://en.wikipedia.org/wiki/FASTA_format
* `-r`, `--vcf-reference <VCF_REFERENCE>` — FASTA file of the sequence the VCF was mapped to (only for vcf input)
* `-t`, `--tree <TREE>` — Name of file containing the tree in newick, nexus, or phylip format.

   If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `--model infer` infers a model from the data. Alternatively, specify the model type. If the specified model requires additional options, use `--model-params` to specify those.

  Default value: `infer`

  Possible values:
  - `infer`:
    Infer GTR parameters from data via Fitch parsimony substitution counts
  - `jc69`
  - `k80`
  - `f81`
  - `hky85`
  - `t92`
  - `tn93`
  - `jtt92`

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as a `key=value` list

   Example: `--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.

   See the exact definitions of the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py
* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`, `joint`

* `--dense <DENSE>` — Use dense representation (stores full probability vectors at each position)

   When combined with `--model infer`, marginal reconstruction runs twice: once to populate profiles for GTR inference, and again with the inferred GTR.

  Possible values: `true`, `false`

* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--zero-based` — Zero-based mutation indexing
* `--reconstruct-tip-states` — Overwrite ambiguous states on tips with the most likely inferred state
* `--report-ambiguous` — Include transitions involving ambiguous states
* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Write augur-compatible node data JSON to this path.

   Contains per-node nucleotide mutations, reconstructed sequences, the alignment mask, genome annotations, and the reference (root) sequence. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.
* `-O`, `--outdir <OUTDIR>` — Directory to write the standard set of output files to
* `--gtr-iterations <GTR_ITERATIONS>` — Number of outer GTR refinement iterations.

   Re-estimates the rate matrix from marginal posterior profiles after each reconstruction pass. Only effective with `--model infer`. Default 0 preserves the current single-pass behavior. Mugration uses 5 by default.

  Default value: `0`
* `--seed <SEED>` — Random seed
* `--sample-from-profile <SAMPLE_FROM_PROFILE>` — How to pick ancestral states from the marginal posterior profile.

   'argmax': most likely state at every node (deterministic, default). 'root': sample from the posterior at the root only, argmax elsewhere (matches augur's `sample_from_profile='root'`). Use `--seed` for reproducible draws. 'all': sample from the posterior at every node.

   Only affects marginal reconstruction (`--method-anc=marginal`).

  Default value: `argmax`

  Possible values: `argmax`, `root`, `all`

* `--const <CONSTANT_SITES>` — Number of constant sites not included in alignment
* `--rescale` — rescale branch lengths
* `--detailed <DETAILED>` — generate a more detailed report
* `--drms <DRMS>` — TSV file containing DRM info. Columns headers: GENOMIC_POSITION, ALT_BASE, DRUG, GENE, SUBSTITUTION
* `-n`, `--num-mut <NUM_MUT>` — number of mutations/nodes that are printed to screen

  Default value: `10`



## `treetime mugration`

Reconstructs discrete ancestral states, for example geographic location, host, or similar. In addition to ancestral states, a GTR model of state transitions is inferred

**Usage:** `treetime mugration [OPTIONS] --attribute <ATTRIBUTE> --metadata <METADATA> --outdir <OUTDIR>`

###### **Options:**

* `-t`, `--tree <TREE>` — Name of file containing the tree in newick, nexus, or phylip format.

   If none is provided, treetime will attempt to build a tree from the alignment using fasttree, iqtree, or raxml (assuming they are installed)
* `--attribute <ATTRIBUTE>` — Attribute to reconstruct, e.g. country
* `-s`, `--metadata <METADATA>` [alias: `states`] — CSV or TSV file with discrete characters. #name,country,continent taxon1,micronesia,oceania ...
* `-w`, `--weights <WEIGHTS>` — CSV or TSV file with probabilities of that a randomly sampled sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.: #country,weight micronesia,0.1 ...
* `--metadata-id-columns <COLUMN>` [alias: `name-column`] — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma and tab.

  Default values: `,`, `	`
* `--output-confidence <OUTPUT_CONFIDENCE>` — Write confidence profile of mugration inference to this path
* `--pc <PC>` — Pseudo-counts. Higher numbers results in 'flatter' models. Default: 1.0
* `--missing-data <MISSING_DATA>` — String indicating missing data

  Default value: `?`
* `--missing-weights-threshold <MISSING_WEIGHTS_THRESHOLD>` — Portion of attribute values that is allowed to not have weights in the weights file

  Default value: `0.5`
* `--iterations <ITERATIONS>` — Number of iterations for GTR model refinement from data

  Default value: `5`
* `--sampling-bias-correction <SAMPLING_BIAS_CORRECTION>` — Rough estimate of how many more events would have been observed if sequences represented an even sample
* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Write augur-compatible node data JSON to this path.

   Contains per-node discrete trait assignments, confidence profiles, entropy, the inferred substitution model, and branch state-change labels. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.
* `-O`, `--outdir <OUTDIR>` — Directory to write the standard set of output files to



## `treetime arg`

Estimates ancestral reassortment graph (ARG)

**Usage:** `treetime arg`



## `treetime schema`

Write JSON Schema definitions for TreeTime data types

**Usage:** `treetime schema [OPTIONS]`

###### **Options:**

* `--for <FOR_FORMAT>` — Which schema to generate

  Default value: `all`

  Possible values: `all`, `version-info`, `progress-event`, `error-response`

* `-o`, `--output <OUTPUT>` — Output file or directory (use "-" for stdout). Directory required when --for=all



## `treetime help-markdown`

Print CLI reference documentation in Markdown format

**Usage:** `treetime help-markdown`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

