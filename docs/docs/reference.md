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
* [`treetime pipeline`↴](#treetime-pipeline)
* [`treetime schema`↴](#treetime-schema)
* [`treetime help-markdown`↴](#treetime-help-markdown)

## `treetime`

Maximum-likelihood phylodynamic inference

Documentation: https://treetime.readthedocs.io/en/stable/
Publication:   https://academic.oup.com/ve/article/4/1/vex042/4794731

**Usage:** `treetime [OPTIONS] <COMMAND>`

Exit status: 0 on success, 1 on failure, 2 on a usage error (invalid flags or config file)

###### **Subcommands:**

* `completions` — Generate shell completions
* `timetree` — Estimates time trees from an initial tree topology, a set of date constraints (e.g. tip dates), and an alignment (optional)
* `optimize` — Optimizes the branch lengths and likelihood of a phylogenetic tree given aligned sequences
* `prune` — Prunes short branches and/or branches without mutations from a phylogenetic tree
* `ancestral` — Reconstructs ancestral sequences and maps mutations to the tree. The `--output-*` flags select the outputs, such as reconstructed sequences, trees, and the inferred GTR model. With marginal reconstruction and `--output-nwk-style=beast` or `nhx`, each branch of the Newick and Nexus trees carries its mutations as a comment, like `A45G,G136T`, with 1-based positions
* `clock` — Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. It will reroot the tree to maximize the clock-like signal and recalculate branch length unless run with `--keep-root`
* `homoplasy` — Reconstructs ancestral sequences and maps mutations to the tree. The tree is then scanned for homoplasies. An excess number of homoplasies might suggest contamination, recombination, culture adaptation or similar
* `mugration` — Reconstructs discrete ancestral states, for example geographic location, host, or similar. In addition to ancestral states, a GTR model of state transitions is inferred
* `pipeline` — Runs an ordered list of analysis commands from one config file, on one dataset, in one process
* `schema` — Write JSON Schema definitions for TreeTime data types
* `help-markdown` — Print CLI reference documentation in Markdown format

###### **Options:**

* `--verbosity <VERBOSITY>` — Set verbosity level of console output

  Default value: `warn`

  Possible values: `off`, `error`, `warn`, `info`, `debug`, `trace`

* `--silent` — Disable all console output. Same as `--verbosity=off`
* `-v`, `--verbose` — Make console output more verbose. Add multiple occurrences to increase verbosity further
* `-q`, `--quiet` — Make console output more quiet. Add multiple occurrences to make output even more quiet
* `--no-progress` — Disable progress bar display
* `--color <WHEN>` — When to color error messages, help, and the log: `auto` colors them when the output is a terminal and the `NO_COLOR` environment variable is not set

  Default value: `auto`

  Possible values: `auto`, `always`, `never`




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

**Usage:** `treetime timetree [OPTIONS]`

###### **Options:**

* `--config <CONFIG>` — Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Aligned FASTA input, plain or compressed (`gz`, `bz2`, `xz`, `zst`); compression is detected by extension.

   Each `--alignment` file is one partition, numbered by its position in the list (1, 2, ...). Use `-` to read uncompressed FASTA from standard input; `-` may appear at most once.
* `-t`, `--tree <TREE>` — Tree in Newick format
* `--tree-dialect <TREE_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `-d`, `--metadata <METADATA>` [alias: `dates`] — Metadata table (CSV or TSV) with one row per sample

   `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the sample, such as its sampling date or its traits.
* `--metadata-id-columns <COLUMN>` [alias: `name-column`] — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma, tab, and semicolon.

  Default values: `,`, `	`, `;`
* `--date-column <COLUMN>` — Label of the column to be used as sampling date (auto-detected when omitted)
* `--date-format <FORMAT>` — Format used to parse string sampling dates (numeric, ISO, and uncertain dates parse regardless)

  Default value: `%Y-%m-%d`
* `--sequence-length <SEQUENCE_LENGTH>` — Length of the sequence, used to calculate expected variation in branch length. Not required if alignment is provided
* `--clock-rate <CLOCK_RATE>` — If specified, the rate of the molecular clock won't be optimized
* `--clock-std-dev <CLOCK_STD_DEV>` — Standard deviation of the provided clock rate estimate
* `--branch-length-mode <BRANCH_LENGTH_MODE>` — How the branch lengths of the time inference are obtained

   `auto` keeps the input branch lengths without an alignment, and also when the longest input branch is longer than 0.1 substitutions per site, because the re-estimated branch lengths are accurate only at short evolutionary distances; otherwise it re-estimates them as `marginal` does. When no input branch has a length, `auto` chooses `marginal`. The run logs the choice.

  Default value: `auto`

  Possible values:
  - `auto`:
    `input` without an alignment or when the longest input branch is longer than 0.1 substitutions per site, `marginal` otherwise
  - `input`:
    Keep the branch lengths of the input tree
  - `marginal`:
    Estimate the branch lengths from the alignment by marginal ancestral reconstruction

* `--confidence` — Estimate the 90% intervals of the node dates

   Runs a final marginal round of the time inference after the iterations and writes the 90% interval of every dated node. The intervals include the clock rate uncertainty when one is known: `--clock-std-dev` for a given rate, or the Gaussian estimate of an estimated rate. The rate uncertainty is added by inferring the dates again at the rate plus and minus one standard deviation and combining both parts in quadrature.
* `--max-grid-points <MAX_GRID_POINTS>` — Largest number of points of one probability grid during time inference.

   A run stops with an error when a grid would need more points. The value bounds single grids, not the total memory of a run. When unset, the run takes the limit of the app or server that runs it, otherwise 1000000.
* `--keep-polytomies` — Keep the polytomies of the input tree

   Without this flag, each polytomy is resolved into binary mergers by the node dates. The resolution is stochastic: the run logs the seed it draws, and `--seed` reproduces the run.
* `--relax <SLACK>` — use an autocorrelated molecular clock. Strength of the gaussian priors on branch specific rate deviation and the coupling of parent and offspring rates can be specified e.g. as --relax 1.0 0.5. Values around 1.0 correspond to weak priors, larger values constrain rate deviations more strongly. Coupling 0 (--relax 1.0 0) corresponds to an un-correlated clock
* `--max-iter <MAX_ITER>` — Maximum number of iterations of the inference loop

   Polytomy resolution and coalescent models of the time inference need at least 2 iterations.

  Default value: `2`
* `--coalescent <TC|opt|skyline>` — Coalescent prior: a fixed time scale in substitutions per site (as in TreeTime 0.x and augur, for example 0.01) or as a time (10y, 26w, 90d); opt fits a constant time scale; skyline fits a piecewise-constant one

   A time is a positive number directly followed by the unit y (years), w (weeks) or d (days). A time scale in substitutions per site is converted to years with the clock rate of each iteration. Without this setting the node times have no coalescent prior.
* `--skyline-n-points <SKYLINE_N_POINTS>` — Number of grid points in skyline coalescent model.

   Only used with --coalescent skyline. Defines how many piecewise linear segments are used to model Tc(t) over time. Must be at least 2. Matches Python v0's default.

  Default value: `20`
* `--skyline-stiffness <SKYLINE_STIFFNESS>` — Smoothing stiffness for the skyline coalescent.

   Penalizes log-fold-changes of the coalescent time scale Tc between adjacent skyline segments: with z = ln Tc, the objective adds `(stiffness/2) * Σ (ln(Tc_{i+1}/Tc_i))^2`. Because it acts on log Tc, the stiffness is dimensionless and scale-independent. Larger values enforce a smoother Tc(t). Only used with --coalescent skyline.

  Default value: `2`
* `--coalescent-confidence <COALESCENT_CONFIDENCE>` — Confidence level for coalescent time scale (Tc) bands, in standard deviations.

   Applies to every inferred coalescent prior (--coalescent opt and --coalescent skyline). The band spans `Tc * exp(±confidence * σ)`, where `σ` is the standard deviation of `ln Tc` from the coalescent likelihood curvature. A fixed time scale is not inferred and therefore has no band.

  Default value: `2`
* `--clock-filter <CLOCK_FILTER>` — Ignore tips that do not follow a loose clock: the number of interquartile ranges from the root-to-tip regression beyond which a tip is an outlier.

   0 switches the filter off.

  Default value: `4`
* `--allow-negative-rate` — Allow negative clock rates.

   By default, rates are forced to be positive. For trees with little temporal signal, removing this restriction gives essentially mid-point rooting.
* `--reroot <REROOT>` — Reroot the tree by temporal-signal optimization.

   Defaults to least-squares when rerooting is enabled. Use --keep-root to keep the input root.

  Possible values: `least-squares`, `min-dev`, `oldest`

* `--reroot-tips <REROOT_TIPS>` — Reroot on the branch leading to a tip, or to the most recent common ancestor of several tips
* `--keep-root` — Keep the root of the input tree instead of rerooting it
* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them.

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

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--partition-rates <PARTITION_RATES>` — Relative substitution rates of the partitions

   Each `--alignment` file is one partition. With two or more partitions, all partitions share the branch lengths, and each partition has one relative rate that scales them for the sites of that partition. `estimate` fits the rate of each partition with the rest of its substitution model fixed (a partition with `--model infer` starts from its fitted rate), then divides all rates by their mean weighted by the number of sites of each partition, so that this mean is 1. The branch lengths are then expected substitutions per site, averaged over all sites of all partitions. `equal` fixes the rate of every partition at 1.

   With one partition, this setting has no effect: the run estimates no relative rate.

  Default value: `estimate`

  Possible values:
  - `estimate`:
    Estimate one relative rate per partition, normalized to a site-weighted mean of 1
  - `equal`:
    Fix the relative rate of every partition at 1

* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

  Default value: `auto`

  Possible values:
  - `auto`:
    Let TreeTime choose the representation; it chooses sparse
  - `dense`:
    Probabilities of every state at every alignment position
  - `sparse`:
    Only the variable positions

* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--include-leaves` — Emit reconstructed leaf (tip) sequences in addition to internal nodes
* `--reconstruct-tip-states` — Same as `--include-leaves --impute-missing-data`: emit tip sequences with ambiguous and unknown states resolved to the most likely inferred state
* `--report-ambiguous` — Include branch mutations from or to the fully ambiguous state (`N` for nucleotides, `X` for amino acids).

   By default these mutations are omitted from the branch mutations of the Newick and Nexus annotations, the Auspice JSON, and the `muts` and `aa_muts` lists of augur node data. Other ambiguity codes, such as `K` or `R`, are always reported. The flag does not affect the MAT outputs, which store `N` as missing data and never contain these mutations.
* `--impute-missing-data` — Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely inferred state.

   Gaps are left as deletions, because they are inferred structure, not missing data. Changes that involve ambiguous characters on terminal branches then disappear from the branch mutations. Only defined for marginal reconstruction; a no-op with a warning under `--method-anc=parsimony`.
* `--no-indels` — Disable indel (insertion and deletion) contributions to branch-length optimization and branch-length distributions.

   When set, branch-length optimization uses the substitution-only likelihood, and the branch distributions of the time inference exclude the Poisson indel term. This matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and TreeTime 0.x. By default, indels count.
* `--divergence-units <DIVERGENCE_UNITS>` — Units of the divergence values in the augur node data and Auspice outputs.

   `mutations-per-site` (default): branch divergence as substitutions per site. `mutations`: count of reconstructed substitutions per branch that change the state: gaps, `N`, and ambiguity codes compatible with the parent state do not count. `mutations` requires reconstructed sequences, so it does not work with input branch lengths (`--branch-length-mode=input`).

  Default value: `mutations-per-site`

  Possible values: `mutations-per-site`, `mutations`

* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Path to output augur-compatible node data JSON.

   Contains per-node dates, branch lengths, clock model parameters, confidence intervals, and divergence metrics. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-gtr <OUTPUT_GTR>` — Path to output GTR model JSON.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-reconstructed-nuc-fasta <OUTPUT_RECONSTRUCTED_NUC_FASTA>` — Path to output reconstructed nucleotide FASTA.

   Holds the reconstructed sequences of the internal nodes, and of the tips with `--include-leaves` or `--reconstruct-tip-states`; `--impute-missing-data` resolves ambiguous tip states. Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).
* `--output-clock-model <OUTPUT_CLOCK_MODEL>` — Path to output clock model JSON.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-clock-csv <OUTPUT_CLOCK_CSV>` — Path to output clock regression CSV.

   One row per sample: the date the regression used, the root-to-tip distance it regressed on, the date the clock model predicts from that distance, and whether the clock filter excluded the sample. Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-confidence-tsv <OUTPUT_CONFIDENCE_TSV>` — Path to output date-confidence-interval TSV.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-tracelog <OUTPUT_TRACELOG>` [alias: `tracelog`] — Path to output iteration-statistics tracelog CSV (monitors convergence).

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-coalescent-tsv <OUTPUT_COALESCENT_TSV>` — Path to output the coalescent time scale as a flat TSV (one row per skyline segment).

   Written when a coalescent prior is set with `--coalescent`. A fixed time scale writes one band-less segment over the tree span. Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-coalescent-csv <OUTPUT_COALESCENT_CSV>` — Path to output the coalescent time scale as a flat CSV (one row per skyline segment).

   Written when a coalescent model is set. Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-coalescent-json <OUTPUT_COALESCENT_JSON>` — Path to output the coalescent time scale as a rich JSON document (inputs + segments).

   Written when a coalescent model is set. Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `-O`, `--output-all <OUTPUT_ALL>` — Write all default output files into this directory.

   Produces the default set of tree and non-tree outputs for the command, using `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which outputs are written.

   Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or supplement the files produced by `--output-all`.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{partition}`, which is replaced by the number. Other tree outputs with `{partition}` in their file name are also written once per file.

   An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats.
* `--output-nwk-style <OUTPUT_NWK_STYLE>` — NWK/Nexus annotation styles to write: `plain`, `beast`, `nhx`.

   Applies to every NWK and Nexus output. With more than one style, files are distinguished by a secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <OUTPUT_TREE_NWK>` — Path to output Newick tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-nexus <OUTPUT_TREE_NEXUS>` — Path to output Nexus tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-auspice <OUTPUT_TREE_AUSPICE>` — Path to output Auspice v2 JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-pb <OUTPUT_TREE_MAT_PB>` — Path to output UShER MAT protobuf tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-json <OUTPUT_TREE_MAT_JSON>` — Path to output UShER MAT JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-graph-json <OUTPUT_TREE_GRAPH_JSON>` — Path to output internal graph JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-dot <OUTPUT_TREE_DOT>` — Path to output Graphviz DOT tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `augur-node-data`, `gtr`, `reconstructed-nuc-fasta`, `clock-model`, `clock-csv`, `confidence-tsv`, `tracelog`, `coalescent-tsv`, `coalescent-csv`, `coalescent-json`

* `--ambiguity-threshold <AMBIGUITY_THRESHOLD>` — Report every state whose posterior probability is at most this far below the probability of the most likely state.

   The value is a difference of probabilities from 0 to 1. With 0, the report names only the states that tie with the most likely state. RAxML-NG uses 0.5 divided by the number of states, for example 0.125 for nucleotides. The report is off when the value is not set.

   The reconstructed FASTA writes a position with more than one reported state as the IUPAC code of these states, or `N` for all four nucleotides. Internal nodes have a report at every position, leaves only at the positions that `--impute-missing-data` imputes.
* `--ambiguity-outputs <AMBIGUITY_OUTPUTS>` — List of the outputs that carry the ambiguity report.

   Without this list, every output of the run that can carry the report carries it. Requires `--ambiguity-threshold`.

  Possible values:
  - `reconstructed-nuc-fasta`:
    Reconstructed nucleotide FASTA

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-target-source <TOPOLOGY_ORDER_TARGET_SOURCE>` — Source for target-order topology sorting

  Possible values: `input`, `reference-topology`, `list`

* `--topology-order-target-file <TOPOLOGY_ORDER_TARGET_FILE>` — File used by list or reference-topology target-order sources
* `--topology-order-target-aggregate <TOPOLOGY_ORDER_TARGET_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `--seed <SEED>` [alias: `rng-seed`] — Random seed

   Without a seed, a run with a random step draws one and logs it, so the run can be reproduced.
* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread
* `--gen-per-year <GEN_PER_YEAR>` — Generations per year for converting the coalescent time scale Tc into an effective population size.

   The coalescent output reports an effective population size `N_e = Tc * gen_per_year`. Tc is already in calendar years, so this factor rescales it into generation units, the standard axis of a skyline plot. Only affects the reported `N_e`; it does not enter the inference.

  Default value: `50`



## `treetime optimize`

Optimizes the branch lengths and likelihood of a phylogenetic tree given aligned sequences

**Usage:** `treetime optimize [OPTIONS]`

###### **Options:**

* `--config <CONFIG>` — Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Aligned FASTA input, plain or compressed (`gz`, `bz2`, `xz`, `zst`); compression is detected by extension.

   Each `--alignment` file is one partition, numbered by its position in the list (1, 2, ...). Use `-` to read uncompressed FASTA from standard input; `-` may appear at most once.
* `-t`, `--tree <TREE>` — Tree in Newick format
* `--tree-dialect <TREE_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them.

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

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--partition-rates <PARTITION_RATES>` — Relative substitution rates of the partitions

   Each `--alignment` file is one partition. With two or more partitions, all partitions share the branch lengths, and each partition has one relative rate that scales them for the sites of that partition. `estimate` fits the rate of each partition with the rest of its substitution model fixed (a partition with `--model infer` starts from its fitted rate), then divides all rates by their mean weighted by the number of sites of each partition, so that this mean is 1. The branch lengths are then expected substitutions per site, averaged over all sites of all partitions. `equal` fixes the rate of every partition at 1.

   With one partition, this setting has no effect: the run estimates no relative rate.

  Default value: `estimate`

  Possible values:
  - `estimate`:
    Estimate one relative rate per partition, normalized to a site-weighted mean of 1
  - `equal`:
    Fix the relative rate of every partition at 1

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

  Default value: `auto`

  Possible values:
  - `auto`:
    Let TreeTime choose the representation; it chooses sparse
  - `dense`:
    Probabilities of every state at every alignment position
  - `sparse`:
    Only the variable positions

* `-O`, `--output-all <OUTPUT_ALL>` — Write all default output files into this directory.

   Produces the default set of tree and non-tree outputs for the command, using `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which outputs are written.

   Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or supplement the files produced by `--output-all`.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{partition}`, which is replaced by the number. Other tree outputs with `{partition}` in their file name are also written once per file.

   An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats.
* `--output-nwk-style <OUTPUT_NWK_STYLE>` — NWK/Nexus annotation styles to write: `plain`, `beast`, `nhx`.

   Applies to every NWK and Nexus output. With more than one style, files are distinguished by a secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <OUTPUT_TREE_NWK>` — Path to output Newick tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-nexus <OUTPUT_TREE_NEXUS>` — Path to output Nexus tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-auspice <OUTPUT_TREE_AUSPICE>` — Path to output Auspice v2 JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-pb <OUTPUT_TREE_MAT_PB>` — Path to output UShER MAT protobuf tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-json <OUTPUT_TREE_MAT_JSON>` — Path to output UShER MAT JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-graph-json <OUTPUT_TREE_GRAPH_JSON>` — Path to output internal graph JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-dot <OUTPUT_TREE_DOT>` — Path to output Graphviz DOT tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--no-indels` — Disable indel (insertion and deletion) contributions to branch-length optimization and branch-length distributions.

   When set, branch-length optimization uses the substitution-only likelihood, and the branch distributions of the time inference exclude the Poisson indel term. This matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and TreeTime 0.x. By default, indels count.
* `--divergence-units <DIVERGENCE_UNITS>` — Units of the divergence values in the augur node data and Auspice outputs.

   `mutations-per-site` (default): branch divergence as substitutions per site. `mutations`: count of reconstructed substitutions per branch that change the state: gaps, `N`, and ambiguity codes compatible with the parent state do not count. `mutations` requires reconstructed sequences, so it does not work with input branch lengths (`--branch-length-mode=input`).

  Default value: `mutations-per-site`

  Possible values: `mutations-per-site`, `mutations`

* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Path to output augur-compatible node data JSON.

   Contains per-node optimized branch lengths (divergence, substitutions per site) and the input alignment and tree paths. The output is compatible with augur export v2 --node-data.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-gtr <OUTPUT_GTR>` — Path to output GTR model JSON.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `augur-node-data`, `gtr`

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-target-source <TOPOLOGY_ORDER_TARGET_SOURCE>` — Source for target-order topology sorting

  Possible values: `input`, `reference-topology`, `list`

* `--topology-order-target-file <TOPOLOGY_ORDER_TARGET_FILE>` — File used by list or reference-topology target-order sources
* `--topology-order-target-aggregate <TOPOLOGY_ORDER_TARGET_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `--max-iter <MAX_ITER>` — Maximum number of iterations of the inference loop

   Polytomy resolution and coalescent models of the time inference need at least 2 iterations.

  Default value: `10`
* `--dp <DP>` — Likelihood convergence threshold. The loop stops when successive likelihoods differ by less than this value, or when a 2-cycle with amplitude below this value is detected

  Default value: `0.1`
* `--damping <DAMPING>` — Damping factor $d$ for outer-loop updates: $b=b_{new}(1-d)+b_{old}d$, where $d=max(damping^{i+1},0.01)$. Higher values reduce oscillation; zero disables damping. Must be in $[0,1)$

  Default value: `0`
* `--branch-length-initial-guess <BRANCH_LENGTH_INITIAL_GUESS>` — Initial branch length estimate before the per-edge optimization.

   - auto: estimate only edges with missing or invalid branch lengths, preserve valid input values (default) - always: estimate all edges, overwriting input branch lengths - never: use input branch lengths as-is; fails if any are missing

  Default value: `auto`

  Possible values:
  - `auto`:
    Estimate only edges with missing or invalid branch lengths, preserve valid input values. No-op when all edges have finite branch lengths
  - `always`:
    Estimate all edges, overwriting input branch lengths
  - `never`:
    Use input branch lengths as-is. Fails if any edge has a missing or invalid branch length

* `--opt-method <OPT_METHOD>` — Per-edge optimizer. Both methods optimize in $\sqrt{t}$ space. Values are `brent` (default) and `newton`

  Default value: `brent`

  Possible values:
  - `brent`:
    Brent's method in $\sqrt{t}$ space
  - `newton`:
    Newton-Raphson in $\sqrt{t}$ space

* `--reroot <REROOT>` — Reroot the tree by minimizing root-to-tip divergence variance.

   By default, optimize keeps the input root. Pass --reroot or --reroot=min-dev to enable divergence-based rerooting. Date-dependent methods (least-squares, oldest) are available in the timetree and clock commands.

  Possible values: `min-dev`

* `--reroot-tips <REROOT_TIPS>` — Reroot on the branch leading to a tip, or to the most recent common ancestor of several tips
* `--keep-root` — Keep the root of the input tree instead of rerooting it
* `--no-collapse-short-branches` — Disable collapsing of internal branches whose optimized length is zero.

   By default the optimize loop contracts internal edges the per-edge optimizer drove to exactly zero that carry no substitutions or indels, turning the resulting binary nodes into polytomies. When set, such edges are kept in the output tree with length zero.
* `--no-merge-siblings` — Disable merging of polytomy siblings that share substitutions.

   By default, sibling branches in a polytomy that carry identical substitutions are grouped under a new internal node. Requires the sparse sequence representation, so this step has no effect under `--representation dense` and the flag is then a no-op.
* `--no-flip-parent-child` — Disable the flip-parent-child step (reversion hoist) in polytomies.

   By default, when a child branch carries the exact reversion of a substitution on the node's parent branch, a new node is inserted that groups the node with that child, removing one mutation per reverted position. When set, those reversions are left in place. Requires the sparse sequence representation, so this step has no effect under `--representation dense` and the flag is then a no-op.
* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime prune`

Prunes short branches and/or branches without mutations from a phylogenetic tree

**Usage:** `treetime prune [OPTIONS]`

###### **Options:**

* `--config <CONFIG>` — Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Aligned FASTA input, plain or compressed (`gz`, `bz2`, `xz`, `zst`); compression is detected by extension.

   Each `--alignment` file is one partition, numbered by its position in the list (1, 2, ...). Use `-` to read uncompressed FASTA from standard input; `-` may appear at most once.
* `-t`, `--tree <TREE>` — Tree in Newick format
* `--tree-dialect <TREE_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `-O`, `--output-all <OUTPUT_ALL>` — Write all default output files into this directory.

   Produces the default set of tree and non-tree outputs for the command, using `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which outputs are written.

   Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or supplement the files produced by `--output-all`.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{partition}`, which is replaced by the number. Other tree outputs with `{partition}` in their file name are also written once per file.

   An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats.
* `--output-nwk-style <OUTPUT_NWK_STYLE>` — NWK/Nexus annotation styles to write: `plain`, `beast`, `nhx`.

   Applies to every NWK and Nexus output. With more than one style, files are distinguished by a secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <OUTPUT_TREE_NWK>` — Path to output Newick tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-nexus <OUTPUT_TREE_NEXUS>` — Path to output Nexus tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-auspice <OUTPUT_TREE_AUSPICE>` — Path to output Auspice v2 JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-pb <OUTPUT_TREE_MAT_PB>` — Path to output UShER MAT protobuf tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-json <OUTPUT_TREE_MAT_JSON>` — Path to output UShER MAT JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-graph-json <OUTPUT_TREE_GRAPH_JSON>` — Path to output internal graph JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-dot <OUTPUT_TREE_DOT>` — Path to output Graphviz DOT tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-gtr <OUTPUT_GTR>` — Path to output GTR model JSON.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `gtr`

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-target-source <TOPOLOGY_ORDER_TARGET_SOURCE>` — Source for target-order topology sorting

  Possible values: `input`, `reference-topology`, `list`

* `--topology-order-target-file <TOPOLOGY_ORDER_TARGET_FILE>` — File used by list or reference-topology target-order sources
* `--topology-order-target-aggregate <TOPOLOGY_ORDER_TARGET_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

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

   Surrounding whitespace is trimmed and empty names are ignored. Use --prune-nodes-list-delimiter to specify a different delimiter.
* `--prune-nodes-list-delimiter <DELIMITER>` — Name separator for `--prune-nodes-list`

   String used to separate node names in the list given to (--prune-nodes-list). Make sure to correctly quote and escape the delimiter according to your shell.

  Default value: `,`
* `-N`, `--prune-nodes-list-file <FILEPATH>` — File containing list of node names to prune

   Path to a file containing node names to remove from the tree, newline-delimited (\n). Surrounding whitespace is trimmed and empty names are ignored.

   Use '-' to read from standard input (stdin).

   Use --prune-nodes-list-file-delimiter to specify a different delimiter.
* `--prune-nodes-list-file-delimiter <DELIMITER>` — Separator for node names in the list file

   Character or string used to separate node names in the list file.

  Default value: `
`
* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime ancestral`

Reconstructs ancestral sequences and maps mutations to the tree. The `--output-*` flags select the outputs, such as reconstructed sequences, trees, and the inferred GTR model. With marginal reconstruction and `--output-nwk-style=beast` or `nhx`, each branch of the Newick and Nexus trees carries its mutations as a comment, like `A45G,G136T`, with 1-based positions

**Usage:** `treetime ancestral [OPTIONS]`

###### **Options:**

* `--config <CONFIG>` — Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Aligned FASTA input, plain or compressed (`gz`, `bz2`, `xz`, `zst`); compression is detected by extension.

   Each `--alignment` file is one partition, numbered by its position in the list (1, 2, ...). Use `-` to read uncompressed FASTA from standard input; `-` may appear at most once.
* `-t`, `--tree <TREE>` — Tree in Newick format
* `--tree-dialect <TREE_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them.

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

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences, which places the mutations on the branches

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

  Default value: `auto`

  Possible values:
  - `auto`:
    Let TreeTime choose the representation; it chooses sparse
  - `dense`:
    Probabilities of every state at every alignment position
  - `sparse`:
    Only the variable positions

* `--sample-from-profile <SAMPLE_FROM_PROFILE>` — How to pick ancestral states from the marginal posterior profile.

   'argmax': most likely state at every node (deterministic, default). 'root': sample from the posterior at the root only, argmax elsewhere (as augur's `sample_from_profile='root'`). 'all': sample from the posterior at every node.

   Only affects marginal reconstruction (`--method-anc=marginal`). Use `--seed` for reproducible draws.

  Default value: `argmax`

  Possible values: `argmax`, `root`, `all`

* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--include-leaves` — Emit reconstructed leaf (tip) sequences in addition to internal nodes
* `--reconstruct-tip-states` — Same as `--include-leaves --impute-missing-data`: emit tip sequences with ambiguous and unknown states resolved to the most likely inferred state
* `--report-ambiguous` — Include branch mutations from or to the fully ambiguous state (`N` for nucleotides, `X` for amino acids).

   By default these mutations are omitted from the branch mutations of the Newick and Nexus annotations, the Auspice JSON, and the `muts` and `aa_muts` lists of augur node data. Other ambiguity codes, such as `K` or `R`, are always reported. The flag does not affect the MAT outputs, which store `N` as missing data and never contain these mutations.
* `--impute-missing-data` — Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely inferred state.

   Gaps are left as deletions, because they are inferred structure, not missing data. Changes that involve ambiguous characters on terminal branches then disappear from the branch mutations. Only defined for marginal reconstruction; a no-op with a warning under `--method-anc=parsimony`.
* `--ignore-missing-alns` — Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.

   Without this flag the run aborts when more than one third of the tips lack a sequence, matching TreeTime v0. Useful when consuming per-CDS translations where some samples have no peptide for a given CDS.
* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Path to output augur-compatible node data JSON.

   Contains per-node nucleotide mutations, reconstructed sequences, the alignment mask, genome annotations, and the reference (root) sequence. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).
* `--output-gtr <OUTPUT_GTR>` — Path to output GTR model JSON.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-reconstructed-nuc-fasta <OUTPUT_RECONSTRUCTED_NUC_FASTA>` — Path to output reconstructed nucleotide FASTA.

   Holds the reconstructed sequences of the internal nodes, and of the tips with `--include-leaves` or `--reconstruct-tip-states`; `--impute-missing-data` resolves ambiguous tip states. Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).
* `--translations <TRANSLATIONS>` — Path template for per-CDS amino-acid FASTA alignments.

   The template must contain a CDS placeholder, replaced with each value from `--cdses` (or each CDS in `--annotation` when `--cdses` is omitted). Both `{cds}` (Nextclade `--output-translations`) and `%GENE` (augur) placeholders are accepted.
* `--cdses <CDS>` [alias: `genes`] — CDS names to reconstruct from `--translations`.

   When omitted, the CDS set is derived from `--annotation`.
* `--annotation <ANNOTATION>` — GFF3 file with CDS coordinates for Augur node data annotations.

   Also supplies the CDS set when `--cdses` is omitted.
* `--aa-root-sequence <AA_ROOT_SEQUENCE>` — FASTA file with one amino-acid root/reference sequence per CDS.

   The root `aa_muts` of augur node data list the differences between this reference and the inferred root, and always omit `X` and gaps.
* `--aa-model <AA_MODEL>` — Amino-acid substitution model. Mirrors the nucleotide `--model`; default `infer` matches augur

  Default value: `infer`

  Possible values:
  - `infer`:
    Infer an amino-acid GTR from the data over the stop-inclusive alphabet. Matches augur
  - `jtt92`:
    Jones-Taylor-Thornton 1992 empirical 20-amino-acid model (no stop codon). Stop codons and any other out-of-alphabet characters in the input are mapped to the unknown state `X`

* `--output-reconstructed-aa-fasta <OUTPUT_RECONSTRUCTED_AA_FASTA>` — Path template for per-CDS reconstructed amino-acid FASTA output (including internal nodes).

   Off by default. When set, the reconstructed sequence of every node is written per CDS. The file name accepts the `{cds}`/`%GENE` placeholders of `--translations`, and needs one when more than one CDS is reconstructed; folder names are used as given.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `-O`, `--output-all <OUTPUT_ALL>` — Write all default output files into this directory.

   Produces the default set of tree and non-tree outputs for the command, using `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which outputs are written.

   Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or supplement the files produced by `--output-all`.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{partition}`, which is replaced by the number. Other tree outputs with `{partition}` in their file name are also written once per file.

   An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats.
* `--output-nwk-style <OUTPUT_NWK_STYLE>` — NWK/Nexus annotation styles to write: `plain`, `beast`, `nhx`.

   Applies to every NWK and Nexus output. With more than one style, files are distinguished by a secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <OUTPUT_TREE_NWK>` — Path to output Newick tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-nexus <OUTPUT_TREE_NEXUS>` — Path to output Nexus tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-auspice <OUTPUT_TREE_AUSPICE>` — Path to output Auspice v2 JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-pb <OUTPUT_TREE_MAT_PB>` — Path to output UShER MAT protobuf tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-json <OUTPUT_TREE_MAT_JSON>` — Path to output UShER MAT JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-graph-json <OUTPUT_TREE_GRAPH_JSON>` — Path to output internal graph JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-dot <OUTPUT_TREE_DOT>` — Path to output Graphviz DOT tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `augur-node-data`, `gtr`, `reconstructed-nuc-fasta`, `reconstructed-aa-fasta`

* `--ambiguity-threshold <AMBIGUITY_THRESHOLD>` — Report every state whose posterior probability is at most this far below the probability of the most likely state.

   The value is a difference of probabilities from 0 to 1. With 0, the report names only the states that tie with the most likely state. RAxML-NG uses 0.5 divided by the number of states, for example 0.125 for nucleotides. The report is off when the value is not set.

   A sequence output writes a position with more than one reported state as the ambiguity code of these states: an IUPAC code for nucleotides, `B`, `Z`, or `J` for amino acids, otherwise `N` or `X`. Internal nodes have a report at every position, leaves only at the positions that `--impute-missing-data` imputes. Nodes whose sequence `--sample-from-profile` draws have no report, and `--method-anc=parsimony` gives no report, because Fitch parsimony computes no posterior.
* `--ambiguity-outputs <AMBIGUITY_OUTPUTS>` — List of the outputs that carry the ambiguity report.

   Without this list, every output of the run that can carry the report carries it. Requires `--ambiguity-threshold`.

  Possible values:
  - `reconstructed-nuc-fasta`:
    Reconstructed nucleotide FASTA
  - `reconstructed-aa-fasta`:
    Reconstructed amino-acid FASTA of each CDS
  - `augur-node-data`:
    Node sequences in the augur node data

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-target-source <TOPOLOGY_ORDER_TARGET_SOURCE>` — Source for target-order topology sorting

  Possible values: `input`, `reference-topology`, `list`

* `--topology-order-target-file <TOPOLOGY_ORDER_TARGET_FILE>` — File used by list or reference-topology target-order sources
* `--topology-order-target-aggregate <TOPOLOGY_ORDER_TARGET_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `--seed <SEED>` [alias: `rng-seed`] — Random seed

   Without a seed, a run with a random step draws one and logs it, so the run can be reproduced.
* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime clock`

Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. It will reroot the tree to maximize the clock-like signal and recalculate branch length unless run with `--keep-root`

**Usage:** `treetime clock [OPTIONS]`

###### **Options:**

* `--config <CONFIG>` — Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-t`, `--tree <TREE>` — Tree in Newick format
* `--tree-dialect <TREE_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `-d`, `--metadata <METADATA>` [alias: `dates`] — Metadata table (CSV or TSV) with one row per sample

   `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the sample, such as its sampling date or its traits.
* `--metadata-id-columns <COLUMN>` [alias: `name-column`] — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma, tab, and semicolon.

  Default values: `,`, `	`, `;`
* `--date-column <COLUMN>` — Label of the column to be used as sampling date (auto-detected when omitted)
* `--date-format <FORMAT>` — Format used to parse string sampling dates (numeric, ISO, and uncertain dates parse regardless)

  Default value: `%Y-%m-%d`
* `--clock-filter <CLOCK_FILTER>` — Ignore tips that do not follow a loose clock: the number of interquartile ranges from the root-to-tip regression beyond which a tip is an outlier.

   0 switches the filter off.

  Default value: `4`
* `--allow-negative-rate` — Allow negative clock rates.

   By default, rates are forced to be positive. For trees with little temporal signal, removing this restriction gives essentially mid-point rooting.
* `--reroot <REROOT>` — Reroot the tree by temporal-signal optimization.

   Defaults to least-squares when rerooting is enabled. Use --keep-root to keep the input root.

  Possible values: `least-squares`, `min-dev`, `oldest`

* `--reroot-tips <REROOT_TIPS>` — Reroot on the branch leading to a tip, or to the most recent common ancestor of several tips
* `--keep-root` — Keep the root of the input tree instead of rerooting it
* `-O`, `--output-all <OUTPUT_ALL>` — Write all default output files into this directory.

   Produces the default set of tree and non-tree outputs for the command, using `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which outputs are written.

   Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or supplement the files produced by `--output-all`.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{partition}`, which is replaced by the number. Other tree outputs with `{partition}` in their file name are also written once per file.

   An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats.
* `--output-nwk-style <OUTPUT_NWK_STYLE>` — NWK/Nexus annotation styles to write: `plain`, `beast`, `nhx`.

   Applies to every NWK and Nexus output. With more than one style, files are distinguished by a secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <OUTPUT_TREE_NWK>` — Path to output Newick tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-nexus <OUTPUT_TREE_NEXUS>` — Path to output Nexus tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-auspice <OUTPUT_TREE_AUSPICE>` — Path to output Auspice v2 JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-pb <OUTPUT_TREE_MAT_PB>` — Path to output UShER MAT protobuf tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-json <OUTPUT_TREE_MAT_JSON>` — Path to output UShER MAT JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-graph-json <OUTPUT_TREE_GRAPH_JSON>` — Path to output internal graph JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-dot <OUTPUT_TREE_DOT>` — Path to output Graphviz DOT tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-clock-model <OUTPUT_CLOCK_MODEL>` — Path to output clock model JSON.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-clock-csv <OUTPUT_CLOCK_CSV>` — Path to output clock regression CSV.

   One row per sample: the date the regression used, the root-to-tip distance it regressed on, the date the clock model predicts from that distance, and whether the clock filter excluded the sample. Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-clock-chart-svg <OUTPUT_CLOCK_CHART_SVG>` — Path to output root-to-tip regression chart in SVG format.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-clock-chart-png <OUTPUT_CLOCK_CHART_PNG>` — Path to output root-to-tip regression chart in PNG format.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-clock-chart-pdf <OUTPUT_CLOCK_CHART_PDF>` — Path to output root-to-tip regression chart in PDF format.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `clock-model`, `clock-csv`, `clock-chart-svg`, `clock-chart-png`, `clock-chart-pdf`

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-target-source <TOPOLOGY_ORDER_TARGET_SOURCE>` — Source for target-order topology sorting

  Possible values: `input`, `reference-topology`, `list`

* `--topology-order-target-file <TOPOLOGY_ORDER_TARGET_FILE>` — File used by list or reference-topology target-order sources
* `--topology-order-target-aggregate <TOPOLOGY_ORDER_TARGET_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime homoplasy`

Reconstructs ancestral sequences and maps mutations to the tree. The tree is then scanned for homoplasies. An excess number of homoplasies might suggest contamination, recombination, culture adaptation or similar

**Usage:** `treetime homoplasy [OPTIONS]`

###### **Options:**

* `--config <CONFIG>` — Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILEPATH>` [alias: `aln`] — Aligned FASTA input, plain or compressed (`gz`, `bz2`, `xz`, `zst`); compression is detected by extension.

   Each `--alignment` file is one partition, numbered by its position in the list (1, 2, ...). Use `-` to read uncompressed FASTA from standard input; `-` may appear at most once.
* `-t`, `--tree <TREE>` — Tree in Newick format
* `--tree-dialect <TREE_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--alphabet <ALPHABET>` — Sequence alphabet

   When omitted, the alphabet is auto-detected from sequence content and falls back to `nuc` when detection is ambiguous.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `-g`, `--model <MODEL>` [alias: `gtr`] — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them.

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

* `--model-params <MODEL_PARAMS>` [alias: `gtr-params`] — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences, which places the mutations on the branches

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

  Default value: `auto`

  Possible values:
  - `auto`:
    Let TreeTime choose the representation; it chooses sparse
  - `dense`:
    Probabilities of every state at every alignment position
  - `sparse`:
    Only the variable positions

* `--gap-fill <GAP_FILL>` — How to handle gap characters in input sequences

   'only-terminal': replace leading and trailing gap characters with the ambiguous character (default, matches v0). 'all': replace all gap characters with the ambiguous character. 'none': leave all gap characters unchanged.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--impute-missing-data` — Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely inferred state.

   Gaps are left as deletions, because they are inferred structure, not missing data. Changes that involve ambiguous characters on terminal branches then disappear from the branch mutations. Only defined for marginal reconstruction; a no-op with a warning under `--method-anc=parsimony`.
* `--ignore-missing-alns` — Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.

   Without this flag the run aborts when more than one third of the tips lack a sequence, matching TreeTime v0. Useful when consuming per-CDS translations where some samples have no peptide for a given CDS.
* `--sample-from-profile <SAMPLE_FROM_PROFILE>` — How to pick ancestral states from the marginal posterior profile.

   'argmax': most likely state at every node (deterministic, default). 'root': sample from the posterior at the root only, argmax elsewhere (as augur's `sample_from_profile='root'`). 'all': sample from the posterior at every node.

   Only affects marginal reconstruction (`--method-anc=marginal`). Use `--seed` for reproducible draws.

  Default value: `argmax`

  Possible values: `argmax`, `root`, `all`

* `--seed <SEED>` [alias: `rng-seed`] — Random seed

   Without a seed, a run with a random step draws one and logs it, so the run can be reproduced.
* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread
* `--const <CONSTANT_SITES>` — Number of constant sites that the alignment leaves out.

   Added to the alignment length to give the number of sites of the genome, which sets the number of sites without mutations and the rate of the Poisson comparison.

  Default value: `0`
* `--rescale <RESCALE>` — Factor that multiplies every branch length of the input tree.

   Use it when the tree is not in substitutions per site, for example `--rescale=0.001` for a tree in substitutions per thousand sites. Scales the reconstruction, the tree outputs, and the reported tree lengths.

  Default value: `1`
* `--detailed` — Add the mutations on terminal branches and the taxa that carry recurrent mutations to the report.

   The statistics JSON always contains them.
* `--drms <DRMS>` — TSV file of drug resistance mutations (DRM) that annotates the report.

   Columns: `GENOMIC_POSITION` (counted from 1), `ALT_BASE`, `DRUG`, `GENE`, `SUBSTITUTION`. A mutation at a listed position gets the gene and the drug, and the substitution when its derived base is a listed `ALT_BASE`.
* `-n`, `--num-mut <NUM_MUT>` — Number of rows of each list in the report.

   The statistics JSON always contains the complete lists.

  Default value: `10`
* `-O`, `--output-all <OUTPUT_ALL>` — Write all default output files into this directory.

   Produces the default set of tree and non-tree outputs for the command, using `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which outputs are written.

   Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or supplement the files produced by `--output-all`.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{partition}`, which is replaced by the number. Other tree outputs with `{partition}` in their file name are also written once per file.

   An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats.
* `--output-nwk-style <OUTPUT_NWK_STYLE>` — NWK/Nexus annotation styles to write: `plain`, `beast`, `nhx`.

   Applies to every NWK and Nexus output. With more than one style, files are distinguished by a secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <OUTPUT_TREE_NWK>` — Path to output Newick tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-nexus <OUTPUT_TREE_NEXUS>` — Path to output Nexus tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-auspice <OUTPUT_TREE_AUSPICE>` — Path to output Auspice v2 JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-pb <OUTPUT_TREE_MAT_PB>` — Path to output UShER MAT protobuf tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-json <OUTPUT_TREE_MAT_JSON>` — Path to output UShER MAT JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-graph-json <OUTPUT_TREE_GRAPH_JSON>` — Path to output internal graph JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-dot <OUTPUT_TREE_DOT>` — Path to output Graphviz DOT tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-homoplasy-stats <OUTPUT_HOMOPLASY_STATS>` — Path to output homoplasy statistics JSON.

   Contains the counts, histograms, Poisson comparison, and complete ranked lists of the report, for substitutions, changes involving ambiguous characters, and insertions and deletions.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-homoplasy-report <OUTPUT_HOMOPLASY_REPORT>` — Path to output homoplasy report text.

   The run also logs the report at info level (`-v`).

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `homoplasy-stats`, `homoplasy-report`

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-target-source <TOPOLOGY_ORDER_TARGET_SOURCE>` — Source for target-order topology sorting

  Possible values: `input`, `reference-topology`, `list`

* `--topology-order-target-file <TOPOLOGY_ORDER_TARGET_FILE>` — File used by list or reference-topology target-order sources
* `--topology-order-target-aggregate <TOPOLOGY_ORDER_TARGET_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`




## `treetime mugration`

Reconstructs discrete ancestral states, for example geographic location, host, or similar. In addition to ancestral states, a GTR model of state transitions is inferred

**Usage:** `treetime mugration [OPTIONS]`

###### **Options:**

* `--config <CONFIG>` — Config file (YAML or JSON) with the settings of the command; `-` reads it from standard input. Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-t`, `--tree <TREE>` — Tree in Newick format
* `--tree-dialect <TREE_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as <structure>,<annotations>. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as x#H1) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast ([&key=value]), nhx ([&&NHX:key=value]) or mrbayes ([&B name value]). For example, enewick,beast reads networks with BEAST annotations and classic,nhx reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--attribute <ATTRIBUTE>` — Attribute to reconstruct, e.g. country
* `-s`, `--metadata <METADATA>` [alias: `states`] — Metadata table (CSV or TSV) with one row per sample

   `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the sample, such as its sampling date or its traits.
* `-w`, `--weights <WEIGHTS>` — CSV or TSV file with probabilities of that a randomly sampled sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.: #country,weight micronesia,0.1 ...
* `--metadata-id-columns <COLUMN>` [alias: `name-column`] — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma, tab, and semicolon.

  Default values: `,`, `	`, `;`
* `--output-confidence-csv <OUTPUT_CONFIDENCE_CSV>` [alias: `confidence`] — Path to output state-probability-profile CSV.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--pc <PC>` — Pseudo-count of the substitution model fit. Higher values give flatter models

  Default value: `1`
* `--missing-data <MISSING_DATA>` — String indicating missing data

  Default value: `?`
* `--missing-weights-threshold <MISSING_WEIGHTS_THRESHOLD>` — Portion of attribute values that is allowed to not have weights in the weights file

  Default value: `0.5`
* `--iterations <ITERATIONS>` — Number of iterations for GTR model refinement from data

  Default value: `5`
* `--sampling-bias-correction <SAMPLING_BIAS_CORRECTION>` — Rough estimate of how many more events would have been observed if sequences represented an even sample
* `--smooth-initial-pi` — Smooth the initial equilibrium frequencies with the pseudo-count before the first reconstruction pass.

   Off by default (TreeTime v0 builds the initial model from raw frequencies and applies the pseudo-count only as infer_gtr regularization). Enabling this flattens the prior for the first pass; it only affects weighted models.
* `--filter-uninformative-root` — Exclude near-uniform root positions from the equilibrium-frequency prior.

   Off by default (TreeTime v0 always folds the root's most-likely state into the prior). Enabling this drops root positions whose posterior carries no phylogenetic signal, removing a state-order-dependent bias at ambiguous roots.
* `--output-augur-node-data <OUTPUT_AUGUR_NODE_DATA>` — Path to output augur-compatible node data JSON.

   Contains per-node discrete trait assignments, confidence profiles, entropy, the inferred substitution model, and branch state-change labels. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-gtr <OUTPUT_GTR>` — Path to output GTR model JSON.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `--output-traits-csv <OUTPUT_TRAITS_CSV>` — Path to output traits CSV.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.
* `-O`, `--output-all <OUTPUT_ALL>` — Write all default output files into this directory.

   Produces the default set of tree and non-tree outputs for the command, using `<dir>/<command>.<ext>` paths. Combine with `--output-selection` to restrict which outputs are written.

   Per-file flags (`--output-tree-nwk`, `--output-augur-node-data`, etc.) override or supplement the files produced by `--output-all`.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, and Newick and Nexus trees in an annotated style. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{partition}`, which is replaced by the number. Other tree outputs with `{partition}` in their file name are also written once per file.

   An output path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats.
* `--output-nwk-style <OUTPUT_NWK_STYLE>` — NWK/Nexus annotation styles to write: `plain`, `beast`, `nhx`.

   Applies to every NWK and Nexus output. With more than one style, files are distinguished by a secondary extension (`.annotated` for beast, `.nhx` for nhx). Default: `plain`.

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <OUTPUT_TREE_NWK>` — Path to output Newick tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-nexus <OUTPUT_TREE_NEXUS>` — Path to output Nexus tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`. With multiple `--output-nwk-style` values, a secondary extension is inserted per style.

   With several `--alignment` files, the styles other than `plain` write one file per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-auspice <OUTPUT_TREE_AUSPICE>` — Path to output Auspice v2 JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-pb <OUTPUT_TREE_MAT_PB>` — Path to output UShER MAT protobuf tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-mat-json <OUTPUT_TREE_MAT_JSON>` — Path to output UShER MAT JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   With several `--alignment` files, one file is written per alignment, and the file name must contain `{partition}`, which is replaced by the number of the alignment (1, 2, ...).

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-graph-json <OUTPUT_TREE_GRAPH_JSON>` — Path to output internal graph JSON tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-tree-dot <OUTPUT_TREE_DOT>` — Path to output Graphviz DOT tree file.

   Takes precedence over paths configured with `--output-all` and `--output-selection`.

   Compression: path ending in `.gz`, `.bz2`, `.xz`, `.zst` writes compressed output. Use `-` to write uncompressed to stdout.

   Parent directories are created if missing.
* `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `augur-node-data`, `gtr`, `confidence-csv`, `traits-csv`

* `--ambiguity-threshold <AMBIGUITY_THRESHOLD>` — Report every state whose posterior probability is at most this far below the probability of the most likely state.

   The value is a difference of probabilities from 0 to 1. With 0, the report names only the states that tie with the most likely state. RAxML-NG uses 0.5 divided by the number of states, for example 0.125 for nucleotides. The report is off when the value is not set.

   A trait output names the reported states of a node with more than one reported state, joined by `|` in the order of the state names, for example `American Samoa|Samoa`. A node with one reported state has an empty cell in the traits CSV and no field in the JSON outputs.
* `--ambiguity-outputs <AMBIGUITY_OUTPUTS>` — List of the outputs that carry the ambiguity report.

   Without this list, every output of the run that can carry the report carries it. Requires `--ambiguity-threshold`.

  Possible values:
  - `traits-csv`:
    Extra column `<attribute>_ambiguous` in the traits CSV
  - `augur-node-data`:
    Extra node field `<attribute>_ambiguous` in the augur node data
  - `auspice`:
    Extra node attribute `<attribute>_ambiguous` in the Auspice JSON

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-target-source <TOPOLOGY_ORDER_TARGET_SOURCE>` — Source for target-order topology sorting

  Possible values: `input`, `reference-topology`, `list`

* `--topology-order-target-file <TOPOLOGY_ORDER_TARGET_FILE>` — File used by list or reference-topology target-order sources
* `--topology-order-target-aggregate <TOPOLOGY_ORDER_TARGET_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime pipeline`

Runs an ordered list of analysis commands from one config file, on one dataset, in one process.

The config file (JSON or YAML) lists named steps, each an analysis command with its arguments. Steps run sequentially; a step may reference an earlier step's outputs with `{{ steps.<name>.outputs.<selection> }}`. Use `--check` to print the resolved plan without running anything, and `--steps` to run only a subset.

**Usage:** `treetime pipeline [OPTIONS] --config <CONFIG>`

###### **Options:**

* `--config <CONFIG>` — Pipeline configuration file (JSON or YAML) describing an ordered list of steps.

   A relative path in the file resolves from the folder of the file, and from the working directory when the file is read from standard input (`-`).
* `--output-all <OUTPUT_ALL>` — Folder of the outputs: each step writes into a subfolder named after the step. Takes precedence over the top-level `output_all` of the pipeline file. A relative path resolves from the working directory
* `--steps <STEPS>` — Run only these named steps, in the pipeline's list order.

   A referenced upstream step that is not selected must already have its outputs on disk.
* `--check` — Resolve and print the plan (steps, inputs, outputs) without running any step
* `-j`, `--jobs <JOBS>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime schema`

Write JSON Schema definitions for TreeTime data types

**Usage:** `treetime schema [OPTIONS]`

###### **Options:**

* `--for <TARGET>` — Which schema to generate

  Default value: `all`

  Possible values: `all`, `version-info`, `progress-event`, `error-response`, `pipeline`, `timetree`, `optimize`, `prune`, `ancestral`, `homoplasy`, `clock`, `mugration`

* `-o`, `--output <OUTPUT>` — Output file or directory (use "-" for stdout). Directory required when --for=all



## `treetime help-markdown`

Print CLI reference documentation in Markdown format

**Usage:** `treetime help-markdown`



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

