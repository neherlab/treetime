---
title: Command-line reference
description: Every command, flag, and default of the treetime command-line program.
---

<!-- WARNING: GENERATED CONTENT -->
<!-- Do not edit the generated file manually! All manual changes will be overwritten by automation. -->
<!-- See developer guide for how to properly work with this. -->

> This page is generated from the help text of the `treetime` program. Report issues at https://github.com/neherlab/treetime/issues




This document contains the automatically generated reference documentation for command-line arguments of the latest version of TreeTime CLI.

If you have TreeTime CLI installed, you can type `treetime --help` to read the latest documentation for your installed version of TreeTime. To generate this document in markdown format, run `treetime help-markdown > reference.md`
  

**Command Overview:**

* [`treetime`↴](#treetime)
* [`treetime timetree`↴](#treetime-timetree)
* [`treetime ancestral`↴](#treetime-ancestral)
* [`treetime clock`↴](#treetime-clock)
* [`treetime mugration`↴](#treetime-mugration)
* [`treetime homoplasy`↴](#treetime-homoplasy)
* [`treetime optimize`↴](#treetime-optimize)
* [`treetime prune`↴](#treetime-prune)
* [`treetime pipeline`↴](#treetime-pipeline)
* [`treetime schema`↴](#treetime-schema)
* [`treetime schema write`↴](#treetime-schema-write)
* [`treetime completions`↴](#treetime-completions)
* [`treetime help-markdown`↴](#treetime-help-markdown)

## `treetime`

Maximum-likelihood phylodynamic inference

Repository:  https://github.com/neherlab/treetime
Publication: https://academic.oup.com/ve/article/4/1/vex042/4794731

**Usage:**

```sh
treetime [OPTIONS] <COMMAND>
```

Exit status: 0 on success, 1 on failure, 2 on a usage error (invalid flags or config file)

###### **Subcommands:**

* `timetree` — Infer a time-scaled tree from a tree, the sampling dates of its tips, and optionally an alignment
* `ancestral` — Reconstruct ancestral sequences and map the mutations to the branches of a tree
* `clock` — Fit a molecular clock by root-to-tip regression, and reroot the tree for the best clock signal
* `mugration` — Reconstruct the ancestral states of a discrete trait, such as location or host, and fit a model of its changes
* `homoplasy` — Reconstruct ancestral sequences and count homoplasies, mutations that occur more than once in a tree
* `optimize` — Optimize the branch lengths of a tree for the likelihood of an alignment
* `prune` — Collapse short branches and branches without mutations, and remove named nodes from a tree
* `pipeline` — Run several analysis commands in order from one config file
* `schema` — Write the JSON schemas of the config files and of the TreeTime data types
* `completions` — Print the shell completions of TreeTime
* `help-markdown` — Print the command-line reference in Markdown

###### **Options:**

* `--verbosity <VERBOSITY>` — Verbosity level of the console output

  Default value: `warn`

  Possible values: `off`, `error`, `warn`, `info`, `debug`, `trace`

* `--silent` — Disable all console output, as `--verbosity=off` does
* `-v`, `--verbose` — Make the console output more verbose; repeat the flag (`-vv`) for more
* `-q`, `--quiet` — Make the console output quieter; repeat the flag (`-qq`) for less
* `--no-progress` — Disable the progress bar
* `--color <WHEN>` — When to color error messages, help, and the log: `auto` colors them when the output is a terminal and the `NO_COLOR` environment variable is not set

  Default value: `auto`

  Possible values: `auto`, `always`, `never`




## `treetime timetree`

Infer a time-scaled tree from a tree, the sampling dates of its tips, and optionally an alignment

**Usage:**

```sh
treetime timetree --tree <FILE> --metadata <FILE> [(--alignment <FILE>)... | --sequence-length <N>] (-O <DIR> | --output-* <FILE>) [OPTIONS]
treetime timetree --config <FILE> [OPTIONS]
```

**Examples:**

```sh
treetime timetree -t tree.nwk -m metadata.tsv -a aln.fasta.xz -O results
treetime timetree -t tree.nwk -m metadata.tsv --sequence-length 10700 --coalescent opt --confidence -O results
treetime timetree --config timetree.yaml -O results
```

###### **Options:**

* `--config <FILE>` — Config file (YAML or JSON) with the settings of the command

   Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILE>` — Aligned FASTA input

   Each `--alignment` file is numbered by its position in the list (1, 2, ...) and is one partition, unless `--partition` entries of `ancestral`, `optimize`, or `timetree` split it into several. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin); `-` may appear at most once.
* `--partition <ENTRY>` — One partition: a set of sites of one `--alignment` file with its own substitution model; repeat the flag for each partition

   Without `--partition`, each `--alignment` file is one partition with the model of `--model`. An entry is a YAML or JSON map with the keys `name` (the name of the partition in outputs and messages), `alignment` (the number of its `--alignment` file, 1, 2, ...; required with two or more files), `sites` (a list of site ranges, numbered from 1: `a`, `a-b`, or `a-b\k` for every k-th site from a to b), `model` (the substitution model, default `infer`) and `model_params` (the `key=value` items of `--model-params` for that model). Each entry carries its own model, so `--model` and `--model-params` cannot be used together with `--partition`. Every alignment needs at least one entry, and every site belongs to exactly one entry. An entry without `sites` takes the sites that the other entries of its alignment leave; one entry per alignment can leave out `sites`. Write site ranges in single quotes or without quotes, because YAML reads a backslash in double quotes as an escape. Example: `--partition "{name: cp1, sites: ['101-1600\3']}" --partition "{name: cp2, sites: ['102-1600\3']}" --partition "{name: cp3, sites: ['103-1600\3'], model: hky85, model_params: ['kappa=4']}" --partition "{name: noncoding}"`.
* `-t`, `--tree <FILE>` — Tree in Newick format

   A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--newick-dialect <NEWICK_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--allow-missing-sequences` — Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.

   Without this flag the run aborts when more than one third of the tips lack a sequence, because the alignment may not belong to the tree. Useful for per-CDS translations, where some samples have no peptide for a given CDS.
* `-m`, `--metadata <FILE>` — Metadata table (CSV or TSV) with one row per sample

   `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the sample, such as its sampling date or its traits. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--metadata-id-columns <COLUMN>` — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma, tab, and semicolon.

  Default values: `,`, `	`, `;`
* `--date-column <COLUMN>` — Label of the column to be used as sampling date (auto-detected when omitted)
* `--date-format <FORMAT>` — Format used to parse string sampling dates (numeric, ISO, and uncertain dates parse regardless)

  Default value: `%Y-%m-%d`
* `--sequence-length <N>` — Length of the sequences, used to calculate the expected variation of the branch lengths

   Required without `--alignment`, and applies only then: with an alignment, the run takes the length from the alignment.
* `--alphabet <ALPHABET>` — Alphabet of the sequences

   `nuc` (default): nucleotides A, C, G and T, with the IUPAC ambiguity codes and `N` for an unknown base. `aa`: the 20 amino acids and the stop codon `*` as a 21st state, with `X` for an unknown residue. `aa-no-stop`: the 20 amino acids without a state for the stop codon.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--clock-rate <CLOCK_RATE>` — Fixed rate of the molecular clock, in substitutions per site per year

   Without this setting, the run estimates the rate.
* `--clock-std-dev <CLOCK_STD_DEV>` — Standard deviation of the provided clock rate estimate

   Applies only with `--confidence`, whose date intervals include the clock rate uncertainty.
* `--relax <SLACK>` — Use an autocorrelated relaxed clock with the given prior strengths, for example `--relax 1.0 0.5`

   SLACK is the strength of the Gaussian prior on the rate deviation of each branch, and COUPLING the strength of the coupling of the rates of parent and child branches. Values around 1.0 give weak priors; larger values constrain the rate deviations more strongly. A COUPLING of 0 (`--relax 1.0 0`) gives an uncorrelated clock.
* `--clock-filter <CLOCK_FILTER>` — Ignore tips that do not follow a loose clock: the number of interquartile ranges from the root-to-tip regression beyond which a tip is an outlier.

   0 switches the filter off.

  Default value: `4`
* `--allow-negative-rate` — Allow negative clock rates.

   By default, rates are forced to be positive. For trees with little temporal signal, removing this restriction gives essentially mid-point rooting.
* `--prune-filtered-tips` — Remove the tips that the clock filter flags from the tree.

   The clock filter flags a tip and leaves it in the tree by default. With this flag, the run removes the flagged tips after the filter: timetree before the time inference, as augur refine does, and clock before the final clock fit. Every output leaves the removed tips out, and a run warning names them. The tips that --keep-tips and --keep-tips-file name stay in the tree, flagged. Requires the clock filter (--clock-filter above 0).
* `--keep-tips <NAME>` — Names of tips that stay in the tree when --prune-filtered-tips removes the flagged tips

   Accepts several names, separated by commas or spaces. A kept tip stays flagged: the clock fit leaves it out, and timetree dates it from the tree. Names that match no tip are ignored.
* `--keep-tips-file <FILE>` — File with names of tips that stay in the tree when --prune-filtered-tips removes the flagged tips, one name per line

   Surrounding whitespace is trimmed and empty lines are ignored. Names that match no tip are ignored, so the file can list samples of several datasets. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--reroot <REROOT>` — Method that chooses the root of the tree

   `least-squares` and `oldest` use the sampling dates; `min-dev` and `midpoint` use only the branch lengths, root once on the input branch lengths, and count every tip. Without a rooting setting the tree is rerooted by least-squares. Use --outgroup to root on named tips, or --keep-root to keep the input root.

  Possible values:
  - `least-squares`:
    Minimizes the residuals of the root-to-tip regression on the sampling dates
  - `min-dev`:
    Minimizes the variance of root-to-tip distances over all tips; more robust than midpoint
  - `oldest`:
    Roots at the middle of the branch above the oldest dated tip
  - `midpoint`:
    Roots at the middle of the longest tip-to-tip path; uses only the two farthest tips

* `--outgroup <OUTGROUP>` — Root the tree on the branch above the most recent common ancestor of the named tips

   With one tip, the root is placed on the branch leading to that tip.
* `--keep-root` — Keep the root of the input tree instead of rerooting it
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

* `--max-iter <N>` — Maximum number of iterations of the inference loop

   Polytomy resolution and coalescent models of the time inference need at least 2 iterations.

  Default value: `2`
* `--confidence` — Estimate the 90% intervals of the node dates

   Runs a final marginal round of the time inference after the iterations and writes the 90% interval of every dated node. The intervals include the clock rate uncertainty when one is known: `--clock-std-dev` for a given rate, or the Gaussian estimate of an estimated rate. The rate uncertainty is added by inferring the dates again at the rate plus and minus one standard deviation and combining both parts in quadrature.
* `--max-grid-points <N>` — Largest number of points of one probability grid during time inference

   A run stops with an error when a grid would need more points. The value bounds single grids, not the total memory of a run. When unset, the run takes the limit of the app or server that runs it, otherwise 1000000.
* `--keep-polytomies` — Keep the polytomies of the input tree

   Without this flag, each polytomy is resolved into binary mergers by the node dates. The resolution is stochastic: the run logs the seed it draws, and `--seed` reproduces the run.
* `--coalescent <TC|opt|skyline>` — Coalescent prior: a fixed time scale in substitutions per site (as in TreeTime 0.x and augur, for example 0.01) or as a time (10y, 26w, 90d); opt fits a constant time scale; skyline fits a piecewise-constant one

   A time is a positive number directly followed by the unit y (years), w (weeks) or d (days). A time scale in substitutions per site is converted to years with the clock rate of each iteration. Without this setting the node times have no coalescent prior.
* `--model-params <MODEL_PARAMS>` — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. In commands with `--partition`, each entry carries its own model and parameters, so `--model-params` is an error next to entries. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--model <MODEL>` — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them. The model applies to every file. In commands with `--partition`, each entry carries its own model, so `--model` is an error next to entries.

   [default: infer]

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

* `--partition-rates <PARTITION_RATES>` — Relative substitution rates of the partitions

   Each `--alignment` file is one partition, or each `--partition` entry when there are entries. With two or more partitions, all partitions share the branch lengths, and each partition has one relative rate that scales them for the sites of that partition. `estimate` fits the rate of each partition with the rest of its substitution model fixed (a partition with `--model infer` starts from its fitted rate), then divides all rates by their mean weighted by the number of sites of each partition, so that this mean is 1. The branch lengths are then expected substitutions per site, averaged over all sites of all partitions. `equal` fixes the rate of every partition at 1.

   The setting applies only with two or more partitions.

   [default: estimate]

  Possible values:
  - `estimate`:
    Estimate one relative rate per partition, normalized to a site-weighted mean of 1
  - `equal`:
    Fix the relative rate of every partition at 1

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

   [default: auto]

  Possible values:
  - `auto`:
    Let TreeTime choose the representation; it chooses sparse
  - `dense`:
    Probabilities of every state at every alignment position
  - `sparse`:
    Only the variable positions

* `--gap-fill <GAP_FILL>` — How to handle the gap characters of the input sequences

   `only-terminal` (default, as in TreeTime 0.x): replace the leading and trailing gaps of each sequence with the unknown character of the alphabet (`N` or `X`). `all`: replace every gap with the unknown character. `none`: keep every gap.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--include-leaves` — Emit reconstructed leaf (tip) sequences in addition to internal nodes
* `--reconstruct-tip-states` — Same as `--include-leaves --impute-missing-data`: emit tip sequences with ambiguous and unknown states resolved to the most likely inferred state
* `--report-ambiguous` — Include branch mutations from or to the fully ambiguous state (`N` for nucleotides, `X` for amino acids).

   By default these mutations are omitted from the branch mutations of the Newick and Nexus annotations, the Auspice JSON, and the `muts` and `aa_muts` lists of augur node data. Other ambiguity codes, such as `K` or `R`, are always reported. The flag does not affect the MAT outputs, which store `N` as missing data and never contain these mutations.
* `--impute-missing-data` — Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely inferred state.

   Gaps are left as deletions, because they are inferred structure, not missing data. Changes that involve ambiguous characters on terminal branches then disappear from the branch mutations. Only defined for marginal reconstruction; a no-op with a warning under `--method-anc=parsimony`.
* `-O`, `--output-all <DIR>` — Folder that receives the default outputs of the command

   Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are created.

   Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs with `{alignment}` in their file name are also written once per file.

   The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}` or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is required when more than one CDS is reconstructed; folder names are used as given.
* `--newick-annotations <NEWICK_ANNOTATIONS>` — Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.

   `plain` writes no comments. `beast` writes a comment such as `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date, and its inferred trait value, as far as the command infers them. Nexus outputs always carry `beast` comments, which FigTree reads.

  Default value: `plain`

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <FILE>` — Path to output Newick tree file.

   `--newick-annotations` sets the node comments.
* `--output-tree-nexus <FILE>` — Path to output Nexus tree file.

   The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
* `--output-tree-auspice <FILE>` — Path to output Auspice v2 JSON tree file
* `--output-tree-mat-pb <FILE>` — Path to output UShER MAT protobuf tree file
* `--output-tree-mat-json <FILE>` — Path to output UShER MAT JSON tree file
* `--output-tree-graph-json <FILE>` — Path to output internal graph JSON tree file
* `--output-tree-dot <FILE>` — Path to output Graphviz DOT tree file
* `-s`, `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `node-data`, `model`, `sequences`, `clock-model`, `rtt`, `dates`, `tracelog`, `coalescent`, `coalescent-json`

* `--no-indels` — Disable indel (insertion and deletion) contributions to branch-length optimization and branch-length distributions.

   When set, branch-length optimization uses the substitution-only likelihood, and the branch distributions of the time inference exclude the Poisson indel term. This matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and TreeTime 0.x. By default, indels count.
* `--divergence-units <DIVERGENCE_UNITS>` — Units of the divergence values in the augur node data and Auspice outputs.

   `mutations-per-site` (default): branch divergence as substitutions per site. `mutations`: count of reconstructed substitutions per branch that change the state: gaps, `N`, and ambiguity codes compatible with the parent state do not count. `mutations` requires reconstructed sequences, so it does not work with input branch lengths (`--branch-length-mode=input`).

  Default value: `mutations-per-site`

  Possible values: `mutations-per-site`, `mutations`

* `--output-node-data <FILE>` — Path to output augur-compatible node data JSON.

   Contains per-node dates, branch lengths, clock model parameters, confidence intervals, and divergence metrics. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.
* `--output-model <FILE>` — Path to output the fitted substitution models as JSON
* `--output-sequences <FILE>` — Path to output the nucleotide sequences of the nodes as FASTA.

   Holds the reconstructed sequences of the internal nodes, and of the tips with `--include-leaves` or `--reconstruct-tip-states`; `--impute-missing-data` resolves ambiguous tip states.
* `--output-clock-model <FILE>` — Path to output clock model JSON
* `--output-rtt <FILE>` — Path to output the root-to-tip regression table.

   One row per sample: the date the regression used, the root-to-tip distance it regressed on, the date the clock model predicts from that distance, and whether the clock filter excluded the sample.
* `--output-dates <FILE>` — Path to output node dates TSV: calendar date and numeric date of every node, with the bounds of the 90% max posterior region with `--confidence`
* `--output-tracelog <FILE>` — Path to output the convergence values of every iteration as a table
* `--output-coalescent <FILE>` — Path to output the coalescent time scale as a table (one row per skyline segment).

   Written when a coalescent prior is set with `--coalescent`. A fixed time scale writes one band-less segment over the tree span.
* `--output-coalescent-json <FILE>` — Path to output the coalescent time scale as a rich JSON document (inputs + segments).

   Written when a coalescent model is set.
* `--ambiguity-threshold <AMBIGUITY_THRESHOLD>` — Report every state whose posterior probability is at most this far below the probability of the most likely state.

   The value is a difference of probabilities from 0 to 1. With 0, the report names only the states that tie with the most likely state. RAxML-NG uses 0.5 divided by the number of states, for example 0.125 for nucleotides. The report is off when the value is not set.

   The reconstructed FASTA writes a position with more than one reported state as the IUPAC code of these states, or `N` for all four nucleotides. Internal nodes have a report at every position, leaves only at the positions that `--impute-missing-data` imputes.
* `--ambiguity-outputs <AMBIGUITY_OUTPUTS>` — List of the outputs that carry the ambiguity report.

   Without this list, every output of the run that can carry the report carries it. Requires `--ambiguity-threshold`.

  Possible values:
  - `sequences`:
    Reconstructed nucleotide FASTA

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-tree <FILE>` — Tree whose leaf order sets the target order

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-list <FILE>` — File with leaf names in the target order, one name per line

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-aggregate <TOPOLOGY_ORDER_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `--seed <SEED>` — Random seed

   Without a seed, a run with a random step draws one and logs it, so the run can be reproduced.
* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread
* `--n-skyline <N>` — Number of grid points of the skyline coalescent model

   Applies only to `--coalescent skyline`. Sets how many piecewise linear segments model Tc(t) over time. Must be at least 2.

   [default: 20]
* `--skyline-stiffness <SKYLINE_STIFFNESS>` — Smoothing stiffness of the skyline coalescent model

   Penalizes log-fold-changes of the coalescent time scale Tc between adjacent skyline segments: with z = ln Tc, the objective adds `(stiffness/2) * Σ (ln(Tc_{i+1}/Tc_i))^2`. Because it acts on log Tc, the stiffness is dimensionless and scale-independent. Larger values enforce a smoother Tc(t). Applies only to `--coalescent skyline`.

   [default: 2.0]
* `--coalescent-confidence <COALESCENT_CONFIDENCE>` — Width of the bands of the coalescent time scale Tc, in standard deviations

   Applies only to the inferred coalescent priors (`--coalescent opt` and `--coalescent skyline`). The band spans `Tc * exp(±confidence * σ)`, where `σ` is the standard deviation of `ln Tc` from the coalescent likelihood curvature. A fixed time scale is not inferred and therefore has no band.

   [default: 2.0]
* `--gen-per-year <GEN_PER_YEAR>` — Generations per year, to convert the coalescent time scale Tc into an effective population size

   The coalescent output reports an effective population size `N_e = Tc * gen_per_year`. Tc is already in calendar years, so this factor rescales it into generation units, the standard axis of a skyline plot. Only affects the reported `N_e`; it does not enter the inference. Applies only with `--coalescent`.

   [default: 50.0]



## `treetime ancestral`

Reconstruct ancestral sequences and map the mutations to the branches of a tree

With marginal reconstruction, each branch of the Nexus trees, and of the Newick trees with `--newick-annotations beast` or `nhx`, carries its mutations as a comment, like `A45G,G136T`, with 1-based positions.

**Usage:**

```sh
treetime ancestral --tree <FILE> (--alignment <FILE>)... (-O <DIR> | --output-* <FILE>) [OPTIONS]
treetime ancestral --config <FILE> [OPTIONS]
```

**Examples:**

```sh
treetime ancestral -t tree.nwk -a aln.fasta.xz -O results
treetime ancestral -t tree.nwk -a ha.fasta.xz -a na.fasta.xz --reconstruct-tip-states --output-sequences 'sequences-{alignment}.fasta'
treetime ancestral -t tree.nwk -a aln.fasta.xz --annotation genome.gff3 --translations 'translations/{cds}.fasta' -O results
```

###### **Options:**

* `--config <FILE>` — Config file (YAML or JSON) with the settings of the command

   Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILE>` — Aligned FASTA input

   Each `--alignment` file is numbered by its position in the list (1, 2, ...) and is one partition, unless `--partition` entries of `ancestral`, `optimize`, or `timetree` split it into several. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin); `-` may appear at most once.
* `--partition <ENTRY>` — One partition: a set of sites of one `--alignment` file with its own substitution model; repeat the flag for each partition

   Without `--partition`, each `--alignment` file is one partition with the model of `--model`. An entry is a YAML or JSON map with the keys `name` (the name of the partition in outputs and messages), `alignment` (the number of its `--alignment` file, 1, 2, ...; required with two or more files), `sites` (a list of site ranges, numbered from 1: `a`, `a-b`, or `a-b\k` for every k-th site from a to b), `model` (the substitution model, default `infer`) and `model_params` (the `key=value` items of `--model-params` for that model). Each entry carries its own model, so `--model` and `--model-params` cannot be used together with `--partition`. Every alignment needs at least one entry, and every site belongs to exactly one entry. An entry without `sites` takes the sites that the other entries of its alignment leave; one entry per alignment can leave out `sites`. Write site ranges in single quotes or without quotes, because YAML reads a backslash in double quotes as an escape. Example: `--partition "{name: cp1, sites: ['101-1600\3']}" --partition "{name: cp2, sites: ['102-1600\3']}" --partition "{name: cp3, sites: ['103-1600\3'], model: hky85, model_params: ['kappa=4']}" --partition "{name: noncoding}"`.
* `-t`, `--tree <FILE>` — Tree in Newick format

   A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--newick-dialect <NEWICK_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--allow-missing-sequences` — Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.

   Without this flag the run aborts when more than one third of the tips lack a sequence, because the alignment may not belong to the tree. Useful for per-CDS translations, where some samples have no peptide for a given CDS.
* `--alphabet <ALPHABET>` — Alphabet of the sequences

   `nuc` (default): nucleotides A, C, G and T, with the IUPAC ambiguity codes and `N` for an unknown base. `aa`: the 20 amino acids and the stop codon `*` as a 21st state, with `X` for an unknown residue. `aa-no-stop`: the 20 amino acids without a state for the stop codon.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--annotation <FILE>` — Genome annotation (GFF3) with the coordinates of the CDSs, for the annotations of the augur node data

   Also supplies the CDSs to reconstruct when `--cds-selection` is not given. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--translations <FILE>` — Path template of the amino-acid alignments, one FASTA file per CDS

   The template must contain a CDS placeholder, which is replaced by each name of `--cds-selection`, or by each CDS of `--annotation` without `--cds-selection`. The placeholder is `{cds}`, as in Nextclade, or `%GENE`, as in augur. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--cds-selection <CDS>` — Names of the CDSs to reconstruct from `--translations`

   Without this setting, the run reconstructs every CDS of `--annotation`.
* `--model-params <MODEL_PARAMS>` — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. In commands with `--partition`, each entry carries its own model and parameters, so `--model-params` is an error next to entries. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--model <MODEL>` — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them. The model applies to every file. In commands with `--partition`, each entry carries its own model, so `--model` is an error next to entries.

   [default: infer]

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

* `--aa-model <AA_MODEL>` — Substitution model of the amino-acid sequences, as `--model` for the nucleotides

   The default `infer`, as in augur, infers the model from the data of each CDS.

  Default value: `infer`

  Possible values:
  - `infer`:
    Infer an amino-acid GTR from the data over the stop-inclusive alphabet. Matches augur
  - `jtt92`:
    Jones-Taylor-Thornton 1992 empirical 20-amino-acid model (no stop codon). Stop codons and any other out-of-alphabet characters in the input are mapped to the unknown state `X`

* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences, which places the mutations on the branches

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

   [default: auto]

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

* `--gap-fill <GAP_FILL>` — How to handle the gap characters of the input sequences

   `only-terminal` (default, as in TreeTime 0.x): replace the leading and trailing gaps of each sequence with the unknown character of the alphabet (`N` or `X`). `all`: replace every gap with the unknown character. `none`: keep every gap.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--aa-root-sequence <FILE>` — FASTA file with one amino-acid reference sequence of the root per CDS

   The root `aa_muts` of augur node data list the differences between this reference and the inferred root, and always omit `X` and gaps. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `-O`, `--output-all <DIR>` — Folder that receives the default outputs of the command

   Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are created.

   Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs with `{alignment}` in their file name are also written once per file.

   The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}` or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is required when more than one CDS is reconstructed; folder names are used as given.
* `--newick-annotations <NEWICK_ANNOTATIONS>` — Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.

   `plain` writes no comments. `beast` writes a comment such as `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date, and its inferred trait value, as far as the command infers them. Nexus outputs always carry `beast` comments, which FigTree reads.

  Default value: `plain`

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <FILE>` — Path to output Newick tree file.

   `--newick-annotations` sets the node comments.
* `--output-tree-nexus <FILE>` — Path to output Nexus tree file.

   The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
* `--output-tree-auspice <FILE>` — Path to output Auspice v2 JSON tree file
* `--output-tree-mat-pb <FILE>` — Path to output UShER MAT protobuf tree file
* `--output-tree-mat-json <FILE>` — Path to output UShER MAT JSON tree file
* `--output-tree-graph-json <FILE>` — Path to output internal graph JSON tree file
* `--output-tree-dot <FILE>` — Path to output Graphviz DOT tree file
* `-s`, `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `node-data`, `model`, `sequences`, `translations`

* `--output-node-data <FILE>` — Path to output augur-compatible node data JSON.

   Contains per-node nucleotide mutations, reconstructed sequences, the alignment mask, genome annotations, and the reference (root) sequence. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.
* `--output-model <FILE>` — Path to output the fitted substitution models as JSON
* `--output-sequences <FILE>` — Path to output the nucleotide sequences of the nodes as FASTA.

   Holds the reconstructed sequences of the internal nodes, and of the tips with `--include-leaves` or `--reconstruct-tip-states`; `--impute-missing-data` resolves ambiguous tip states.
* `--output-translations <FILE>` — Path to output the amino-acid sequences of the nodes as FASTA, one file per CDS

   Holds the reconstructed sequence of every node, internal nodes included.
* `--ambiguity-threshold <AMBIGUITY_THRESHOLD>` — Report every state whose posterior probability is at most this far below the probability of the most likely state.

   The value is a difference of probabilities from 0 to 1. With 0, the report names only the states that tie with the most likely state. RAxML-NG uses 0.5 divided by the number of states, for example 0.125 for nucleotides. The report is off when the value is not set.

   A sequence output writes a position with more than one reported state as the ambiguity code of these states: an IUPAC code for nucleotides, `B`, `Z`, or `J` for amino acids, otherwise `N` or `X`. Internal nodes have a report at every position, leaves only at the positions that `--impute-missing-data` imputes. Nodes whose sequence `--sample-from-profile` draws have no report. The setting applies only when a marginal reconstruction runs (`--method-anc marginal`, or the amino-acid reconstruction of `--translations`), because Fitch parsimony computes no posterior.
* `--ambiguity-outputs <AMBIGUITY_OUTPUTS>` — List of the outputs that carry the ambiguity report.

   Without this list, every output of the run that can carry the report carries it. Requires `--ambiguity-threshold`.

  Possible values:
  - `sequences`:
    Reconstructed nucleotide FASTA
  - `translations`:
    Reconstructed amino-acid FASTA of each CDS
  - `node-data`:
    Node sequences in the augur node data

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-tree <FILE>` — Tree whose leaf order sets the target order

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-list <FILE>` — File with leaf names in the target order, one name per line

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-aggregate <TOPOLOGY_ORDER_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `--seed <SEED>` — Random seed

   Without a seed, a run with a random step draws one and logs it, so the run can be reproduced.
* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread
* `--include-leaves` — Emit reconstructed leaf (tip) sequences in addition to internal nodes
* `--impute-missing-data` — Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely inferred state.

   Gaps are left as deletions, because they are inferred structure, not missing data. Changes that involve ambiguous characters on terminal branches then disappear from the branch mutations. Only defined for marginal reconstruction; a no-op with a warning under `--method-anc=parsimony`.
* `--reconstruct-tip-states` — Same as `--include-leaves --impute-missing-data`: emit tip sequences with ambiguous and unknown states resolved to the most likely inferred state
* `--report-ambiguous` — Include branch mutations from or to the fully ambiguous state (`N` for nucleotides, `X` for amino acids).

   By default these mutations are omitted from the branch mutations of the Newick and Nexus annotations, the Auspice JSON, and the `muts` and `aa_muts` lists of augur node data. Other ambiguity codes, such as `K` or `R`, are always reported. The flag does not affect the MAT outputs, which store `N` as missing data and never contain these mutations.



## `treetime clock`

Fit a molecular clock by root-to-tip regression, and reroot the tree for the best clock signal

The regression quantifies how clock-like the tree is. Unless `--keep-root` is given, the run reroots the tree to maximize the clock signal.

**Usage:**

```sh
treetime clock --tree <FILE> --metadata <FILE> (-O <DIR> | --output-* <FILE>) [OPTIONS]
treetime clock --config <FILE> [OPTIONS]
```

**Examples:**

```sh
treetime clock -t tree.nwk -m metadata.tsv -O results
treetime clock -t tree.nwk -m metadata.tsv --keep-root --output-rtt rtt.tsv
```

###### **Options:**

* `--config <FILE>` — Config file (YAML or JSON) with the settings of the command

   Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-t`, `--tree <FILE>` — Tree in Newick format

   A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--newick-dialect <NEWICK_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `-m`, `--metadata <FILE>` — Metadata table (CSV or TSV) with one row per sample

   `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the sample, such as its sampling date or its traits. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--metadata-id-columns <COLUMN>` — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

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
* `--prune-filtered-tips` — Remove the tips that the clock filter flags from the tree.

   The clock filter flags a tip and leaves it in the tree by default. With this flag, the run removes the flagged tips after the filter: timetree before the time inference, as augur refine does, and clock before the final clock fit. Every output leaves the removed tips out, and a run warning names them. The tips that --keep-tips and --keep-tips-file name stay in the tree, flagged. Requires the clock filter (--clock-filter above 0).
* `--keep-tips <NAME>` — Names of tips that stay in the tree when --prune-filtered-tips removes the flagged tips

   Accepts several names, separated by commas or spaces. A kept tip stays flagged: the clock fit leaves it out, and timetree dates it from the tree. Names that match no tip are ignored.
* `--keep-tips-file <FILE>` — File with names of tips that stay in the tree when --prune-filtered-tips removes the flagged tips, one name per line

   Surrounding whitespace is trimmed and empty lines are ignored. Names that match no tip are ignored, so the file can list samples of several datasets. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--reroot <REROOT>` — Method that chooses the root of the tree

   `least-squares` and `oldest` use the sampling dates; `min-dev` and `midpoint` use only the branch lengths, root once on the input branch lengths, and count every tip. Without a rooting setting the tree is rerooted by least-squares. Use --outgroup to root on named tips, or --keep-root to keep the input root.

  Possible values:
  - `least-squares`:
    Minimizes the residuals of the root-to-tip regression on the sampling dates
  - `min-dev`:
    Minimizes the variance of root-to-tip distances over all tips; more robust than midpoint
  - `oldest`:
    Roots at the middle of the branch above the oldest dated tip
  - `midpoint`:
    Roots at the middle of the longest tip-to-tip path; uses only the two farthest tips

* `--outgroup <OUTGROUP>` — Root the tree on the branch above the most recent common ancestor of the named tips

   With one tip, the root is placed on the branch leading to that tip.
* `--keep-root` — Keep the root of the input tree instead of rerooting it
* `-O`, `--output-all <DIR>` — Folder that receives the default outputs of the command

   Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are created.

   Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs with `{alignment}` in their file name are also written once per file.

   The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}` or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is required when more than one CDS is reconstructed; folder names are used as given.
* `--newick-annotations <NEWICK_ANNOTATIONS>` — Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.

   `plain` writes no comments. `beast` writes a comment such as `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date, and its inferred trait value, as far as the command infers them. Nexus outputs always carry `beast` comments, which FigTree reads.

  Default value: `plain`

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <FILE>` — Path to output Newick tree file.

   `--newick-annotations` sets the node comments.
* `--output-tree-nexus <FILE>` — Path to output Nexus tree file.

   The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
* `--output-tree-auspice <FILE>` — Path to output Auspice v2 JSON tree file
* `--output-tree-mat-pb <FILE>` — Path to output UShER MAT protobuf tree file
* `--output-tree-mat-json <FILE>` — Path to output UShER MAT JSON tree file
* `--output-tree-graph-json <FILE>` — Path to output internal graph JSON tree file
* `--output-tree-dot <FILE>` — Path to output Graphviz DOT tree file
* `-s`, `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `clock-model`, `rtt`, `rtt-chart`

* `--output-clock-model <FILE>` — Path to output clock model JSON
* `--output-rtt <FILE>` — Path to output the root-to-tip regression table.

   One row per sample: the date the regression used, the root-to-tip distance it regressed on, the date the clock model predicts from that distance, and whether the clock filter excluded the sample.
* `--output-rtt-chart <FILE>` — Path to output the root-to-tip regression chart.

   The chart plots the root-to-tip divergence of each dated sample against its date, with the least-squares fit, its extrapolation to the root, and the samples that the clock filter excluded. It states the clock rate, the root date and R² of the fit. The extension chooses the image format: `.svg`, `.png` or `.pdf`.
* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-tree <FILE>` — Tree whose leaf order sets the target order

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-list <FILE>` — File with leaf names in the target order, one name per line

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-aggregate <TOPOLOGY_ORDER_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime mugration`

Reconstruct the ancestral states of a discrete trait, such as location or host, and fit a model of its changes

**Usage:**

```sh
treetime mugration --tree <FILE> --metadata <FILE> --attribute <COLUMN> (-O <DIR> | --output-* <FILE>) [OPTIONS]
treetime mugration --config <FILE> [OPTIONS]
```

**Examples:**

```sh
treetime mugration -t tree.nwk -m metadata.tsv --attribute country -O results
treetime mugration -t tree.nwk -m metadata.tsv --attribute country --weights weights.tsv --confidence -O results
```

###### **Options:**

* `--config <FILE>` — Config file (YAML or JSON) with the settings of the command

   Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-t`, `--tree <FILE>` — Tree in Newick format

   A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--newick-dialect <NEWICK_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `-m`, `--metadata <FILE>` — Metadata table (CSV or TSV) with one row per sample

   `--metadata-id-columns` names the column that links a row to a tree tip. The other columns hold the data of the sample, such as its sampling date or its traits. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--metadata-id-columns <COLUMN>` — Candidate column name(s) holding the taxon identifier that links metadata to tree tips

   The first listed column that is present in the header is used. Matching is case-insensitive.

  Default values: `strain`, `name`, `accession`
* `--metadata-delimiters <CHAR>` — Candidate field delimiter(s) for the metadata table

   The delimiter actually present in the file is used. Defaults to comma, tab, and semicolon.

  Default values: `,`, `	`, `;`
* `--attribute <COLUMN>` — Metadata column with the discrete trait to reconstruct, for example `country`
* `--weights <FILE>` — Table (CSV or TSV) with the equilibrium probability of each state

   The weight of a state is the probability that a randomly sampled sequence at equilibrium has the state, for example the population of each country. The first column holds the state and the second its weight, such as the row `micronesia,0.1` under the header `country,weight`. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--lat-longs <FILE>` — Table with the coordinates of the trait values, which give the Auspice tree a map

   The format of augur's `--lat-longs`: one location per line, `<resolution> <name> <latitude> <longitude>`, with the fields separated by tabs, or by white space in a line without a tab. `<resolution>` is the `--attribute` name, and `<name>` a value of the trait; both match in lower case. A later line for the same location replaces an earlier one, and lines that start with `#` are skipped. Latitudes are in -90 to 90 degrees, longitudes in -180 to 180 degrees. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--missing-data <MISSING_DATA>` — Value of the trait column that marks missing data

  Default value: `?`
* `--pc <PC>` — Pseudo-count of the substitution model fit. Higher values give flatter models

  Default value: `1`
* `--missing-weights-threshold <MISSING_WEIGHTS_THRESHOLD>` — Largest fraction of the trait values that may have no weight in the `--weights` file

  Default value: `0.5`
* `--model-iterations <N>` — Number of iterations of the refinement of the substitution model from the data

  Default value: `5`
* `--sampling-bias-correction <SAMPLING_BIAS_CORRECTION>` — Correction for biased sampling: a rough estimate of how many more state changes the tree would show if the samples were an even sample of the population

   The run multiplies the fitted rate of the substitution model by this factor, for example 2.5 for an estimate of 2.5 times more changes. Without this setting, the run applies no correction.
* `--smooth-initial-pi` — Smooth the initial equilibrium frequencies with the pseudo-count before the first reconstruction pass

   By default, the initial model takes the frequencies of `--weights` as they are, and the pseudo-count applies only to the fits of the model. This setting flattens the prior of the first pass. It affects only models with `--weights`.
* `--filter-uninformative-root` — Exclude a near-uniform root posterior from the prior of the equilibrium frequencies

   By default, the most likely state of the root always enters the prior. This setting leaves out a root whose posterior carries no phylogenetic signal, which removes a bias at ambiguous roots that depends on the order of the states.
* `-O`, `--output-all <DIR>` — Folder that receives the default outputs of the command

   Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are created.

   Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs with `{alignment}` in their file name are also written once per file.

   The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}` or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is required when more than one CDS is reconstructed; folder names are used as given.
* `--newick-annotations <NEWICK_ANNOTATIONS>` — Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.

   `plain` writes no comments. `beast` writes a comment such as `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date, and its inferred trait value, as far as the command infers them. Nexus outputs always carry `beast` comments, which FigTree reads.

  Default value: `plain`

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <FILE>` — Path to output Newick tree file.

   `--newick-annotations` sets the node comments.
* `--output-tree-nexus <FILE>` — Path to output Nexus tree file.

   The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
* `--output-tree-auspice <FILE>` — Path to output Auspice v2 JSON tree file
* `--output-tree-mat-pb <FILE>` — Path to output UShER MAT protobuf tree file
* `--output-tree-mat-json <FILE>` — Path to output UShER MAT JSON tree file
* `--output-tree-graph-json <FILE>` — Path to output internal graph JSON tree file
* `--output-tree-dot <FILE>` — Path to output Graphviz DOT tree file
* `-s`, `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `node-data`, `model`, `trait-probabilities`, `traits`

* `--confidence` — Write the probability of every state at every node

   Adds the `trait-probabilities` table to the outputs of `--output-all`, also when `--output-selection` lists other outputs. Without `--output-all`, give the file with `--output-trait-probabilities`.
* `--output-node-data <FILE>` — Path to output augur-compatible node data JSON.

   Contains per-node discrete trait assignments, confidence profiles, entropy, the inferred substitution model, and branch state-change labels. The output is compatible with augur export v2 --node-data for Nextstrain pipeline integration.
* `--output-model <FILE>` — Path to output the fitted substitution models as JSON
* `--output-traits <FILE>` — Path to output the inferred trait state of every node as a table
* `--output-trait-probabilities <FILE>` — Path to output the state probabilities of every node as a table
* `--ambiguity-threshold <AMBIGUITY_THRESHOLD>` — Report every state whose posterior probability is at most this far below the probability of the most likely state.

   The value is a difference of probabilities from 0 to 1. With 0, the report names only the states that tie with the most likely state. RAxML-NG uses 0.5 divided by the number of states, for example 0.125 for nucleotides. The report is off when the value is not set.

   A trait output names the reported states of a node with more than one reported state, joined by `|` in the order of the state names, for example `American Samoa|Samoa`. A node with one reported state has an empty cell in the traits table and no field in the JSON outputs.
* `--ambiguity-outputs <AMBIGUITY_OUTPUTS>` — List of the outputs that carry the ambiguity report.

   Without this list, every output of the run that can carry the report carries it. Requires `--ambiguity-threshold`.

  Possible values:
  - `traits`:
    Extra column `<attribute>_ambiguous` in the traits table
  - `node-data`:
    Extra node field `<attribute>_ambiguous` in the augur node data
  - `auspice`:
    Extra node attribute `<attribute>_ambiguous` in the Auspice JSON

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-tree <FILE>` — Tree whose leaf order sets the target order

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-list <FILE>` — File with leaf names in the target order, one name per line

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-aggregate <TOPOLOGY_ORDER_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime homoplasy`

Reconstruct ancestral sequences and count homoplasies, mutations that occur more than once in a tree

An excess of homoplasies can point to contamination, recombination, culture adaptation, or similar processes.

**Usage:**

```sh
treetime homoplasy --tree <FILE> --alignment <FILE> (-O <DIR> | --output-* <FILE>) [OPTIONS]
treetime homoplasy --config <FILE> [OPTIONS]
```

**Examples:**

```sh
treetime homoplasy -t tree.nwk -a aln.fasta.xz -O results
treetime homoplasy -t tree.nwk -a aln.fasta.xz --drms drms.tsv --detailed --output-homoplasy-report report.txt
```

###### **Options:**

* `--config <FILE>` — Config file (YAML or JSON) with the settings of the command

   Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILE>` — Aligned FASTA input

   Each `--alignment` file is numbered by its position in the list (1, 2, ...) and is one partition, unless `--partition` entries of `ancestral`, `optimize`, or `timetree` split it into several. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin); `-` may appear at most once.
* `-t`, `--tree <FILE>` — Tree in Newick format

   A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--newick-dialect <NEWICK_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--allow-missing-sequences` — Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.

   Without this flag the run aborts when more than one third of the tips lack a sequence, because the alignment may not belong to the tree. Useful for per-CDS translations, where some samples have no peptide for a given CDS.
* `--alphabet <ALPHABET>` — Alphabet of the sequences

   `nuc` (default): nucleotides A, C, G and T, with the IUPAC ambiguity codes and `N` for an unknown base. `aa`: the 20 amino acids and the stop codon `*` as a 21st state, with `X` for an unknown residue. `aa-no-stop`: the 20 amino acids without a state for the stop codon.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--model-params <MODEL_PARAMS>` — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. In commands with `--partition`, each entry carries its own model and parameters, so `--model-params` is an error next to entries. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--model <MODEL>` — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them. The model applies to every file. In commands with `--partition`, each entry carries its own model, so `--model` is an error next to entries.

   [default: infer]

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

* `--method-anc <METHOD_ANC>` — Method used for reconstructing ancestral sequences, which places the mutations on the branches

  Default value: `marginal`

  Possible values: `marginal`, `parsimony`

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

   [default: auto]

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

* `--gap-fill <GAP_FILL>` — How to handle the gap characters of the input sequences

   `only-terminal` (default, as in TreeTime 0.x): replace the leading and trailing gaps of each sequence with the unknown character of the alphabet (`N` or `X`). `all`: replace every gap with the unknown character. `none`: keep every gap.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--impute-missing-data` — Resolve ambiguous and unknown tip states (`N` and IUPAC codes such as `R`) to the most likely inferred state.

   Gaps are left as deletions, because they are inferred structure, not missing data. Changes that involve ambiguous characters on terminal branches then disappear from the branch mutations. Only defined for marginal reconstruction; a no-op with a warning under `--method-anc=parsimony`.
* `--const <N>` — Number of constant sites that the alignment leaves out.

   Added to the alignment length to give the number of sites of the genome, which sets the number of sites without mutations and the rate of the Poisson comparison.

  Default value: `0`
* `--rescale <RESCALE>` — Factor that multiplies every branch length of the input tree.

   Use it when the tree is not in substitutions per site, for example `--rescale=0.001` for a tree in substitutions per thousand sites. Scales the reconstruction, the tree outputs, and the reported tree lengths.

  Default value: `1`
* `--detailed` — Add the mutations on terminal branches and the taxa that carry recurrent mutations to the report.

   The statistics JSON always contains them.
* `--drms <FILE>` — TSV file of drug resistance mutations (DRM) that annotates the report.

   Columns: `GENOMIC_POSITION` (counted from 1), `ALT_BASE`, `DRUG`, `GENE`, `SUBSTITUTION`. A mutation at a listed position gets the gene and the drug, and the substitution when its derived base is a listed `ALT_BASE`. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `-n`, `--num-mut <N>` — Number of rows of each list in the report.

   The statistics JSON always contains the complete lists.

  Default value: `10`
* `-O`, `--output-all <DIR>` — Folder that receives the default outputs of the command

   Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are created.

   Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs with `{alignment}` in their file name are also written once per file.

   The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}` or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is required when more than one CDS is reconstructed; folder names are used as given.
* `--newick-annotations <NEWICK_ANNOTATIONS>` — Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.

   `plain` writes no comments. `beast` writes a comment such as `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date, and its inferred trait value, as far as the command infers them. Nexus outputs always carry `beast` comments, which FigTree reads.

  Default value: `plain`

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <FILE>` — Path to output Newick tree file.

   `--newick-annotations` sets the node comments.
* `--output-tree-nexus <FILE>` — Path to output Nexus tree file.

   The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
* `--output-tree-auspice <FILE>` — Path to output Auspice v2 JSON tree file
* `--output-tree-mat-pb <FILE>` — Path to output UShER MAT protobuf tree file
* `--output-tree-mat-json <FILE>` — Path to output UShER MAT JSON tree file
* `--output-tree-graph-json <FILE>` — Path to output internal graph JSON tree file
* `--output-tree-dot <FILE>` — Path to output Graphviz DOT tree file
* `-s`, `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `homoplasy-stats`, `homoplasy-report`

* `--output-homoplasy-stats <FILE>` — Path to output homoplasy statistics JSON.

   Contains the counts, histograms, Poisson comparison, and complete ranked lists of the report, for substitutions, changes involving ambiguous characters, and insertions and deletions.
* `--output-homoplasy-report <FILE>` — Path to output homoplasy report text.

   The run also logs the report at info level (`-v`).
* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-tree <FILE>` — Tree whose leaf order sets the target order

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-list <FILE>` — File with leaf names in the target order, one name per line

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-aggregate <TOPOLOGY_ORDER_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `--seed <SEED>` — Random seed

   Without a seed, a run with a random step draws one and logs it, so the run can be reproduced.
* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime optimize`

Optimize the branch lengths of a tree for the likelihood of an alignment

**Usage:**

```sh
treetime optimize --tree <FILE> (--alignment <FILE>)... (-O <DIR> | --output-* <FILE>) [OPTIONS]
treetime optimize --config <FILE> [OPTIONS]
```

**Examples:**

```sh
treetime optimize -t tree.nwk -a aln.fasta.xz -O results
treetime optimize -t tree.nwk -a aln.fasta.xz --model hky85 --reroot min-dev --output-tree-nwk optimized.nwk
```

###### **Options:**

* `--config <FILE>` — Config file (YAML or JSON) with the settings of the command

   Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILE>` — Aligned FASTA input

   Each `--alignment` file is numbered by its position in the list (1, 2, ...) and is one partition, unless `--partition` entries of `ancestral`, `optimize`, or `timetree` split it into several. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin); `-` may appear at most once.
* `--partition <ENTRY>` — One partition: a set of sites of one `--alignment` file with its own substitution model; repeat the flag for each partition

   Without `--partition`, each `--alignment` file is one partition with the model of `--model`. An entry is a YAML or JSON map with the keys `name` (the name of the partition in outputs and messages), `alignment` (the number of its `--alignment` file, 1, 2, ...; required with two or more files), `sites` (a list of site ranges, numbered from 1: `a`, `a-b`, or `a-b\k` for every k-th site from a to b), `model` (the substitution model, default `infer`) and `model_params` (the `key=value` items of `--model-params` for that model). Each entry carries its own model, so `--model` and `--model-params` cannot be used together with `--partition`. Every alignment needs at least one entry, and every site belongs to exactly one entry. An entry without `sites` takes the sites that the other entries of its alignment leave; one entry per alignment can leave out `sites`. Write site ranges in single quotes or without quotes, because YAML reads a backslash in double quotes as an escape. Example: `--partition "{name: cp1, sites: ['101-1600\3']}" --partition "{name: cp2, sites: ['102-1600\3']}" --partition "{name: cp3, sites: ['103-1600\3'], model: hky85, model_params: ['kappa=4']}" --partition "{name: noncoding}"`.
* `-t`, `--tree <FILE>` — Tree in Newick format

   A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--newick-dialect <NEWICK_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--allow-missing-sequences` — Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.

   Without this flag the run aborts when more than one third of the tips lack a sequence, because the alignment may not belong to the tree. Useful for per-CDS translations, where some samples have no peptide for a given CDS.
* `--alphabet <ALPHABET>` — Alphabet of the sequences

   `nuc` (default): nucleotides A, C, G and T, with the IUPAC ambiguity codes and `N` for an unknown base. `aa`: the 20 amino acids and the stop codon `*` as a 21st state, with `X` for an unknown residue. `aa-no-stop`: the 20 amino acids without a state for the stop codon.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--reroot <REROOT>` — Method that chooses the root of the tree from the branch lengths

   The tree is rerooted once, on the branch lengths after the initial guess, and every tip counts. Without this setting, the input root stays. The rerooting methods that need dates (`least-squares`, `oldest`) are available where the inputs have dates.

  Possible values:
  - `min-dev`:
    Minimizes the variance of root-to-tip distances over all tips; more robust than midpoint
  - `midpoint`:
    Roots at the middle of the longest tip-to-tip path; uses only the two farthest tips

* `--outgroup <OUTGROUP>` — Root the tree on the branch above the most recent common ancestor of the named tips

   With one tip, the root is placed on the branch leading to that tip.
* `--keep-root` — Keep the root of the input tree instead of rerooting it
* `--max-iter <N>` — Maximum number of iterations of the inference loop

   Polytomy resolution and coalescent models of the time inference need at least 2 iterations.

  Default value: `10`
* `--tolerance <TOLERANCE>` — Convergence threshold of the log-likelihood

   The loop stops when successive log-likelihoods differ by less than this value, or when it detects a 2-cycle with an amplitude below this value.

  Default value: `0.1`
* `--damping <DAMPING>` — Damping of the branch length updates of the loop, from 0 (inclusive) to 1 (exclusive)

   Iteration i updates each branch length b to `b = b_new * (1 - d) + b_old * d`, with `d = max(damping^(i+1), 0.01)`. Higher values reduce oscillation; 0 switches damping off.

  Default value: `0`
* `--initial-branch-lengths <INITIAL_BRANCH_LENGTHS>` — Initial branch lengths before the optimization of each branch

  Default value: `auto`

  Possible values:
  - `auto`:
    Estimate only the branches with a missing or invalid length, and keep the valid input lengths
  - `always`:
    Estimate every branch length, replacing the input lengths
  - `never`:
    Keep the input branch lengths; the run fails when a branch has a missing or invalid length

* `--optimizer <OPTIMIZER>` — Optimizer of each branch length; both methods optimize sqrt(t) in place of the branch length t

  Default value: `brent`

  Possible values:
  - `brent`:
    Brent's method on sqrt(t)
  - `newton`:
    Newton-Raphson on sqrt(t)

* `--no-collapse-zero-branches` — Disable collapsing of internal branches whose optimized length is zero.

   By default the optimize loop contracts internal edges the per-edge optimizer drove to exactly zero that carry no substitutions or indels, turning the resulting binary nodes into polytomies. When set, such edges are kept in the output tree with length zero.
* `--no-merge-shared-mutations` — Disable merging of polytomy siblings that share substitutions.

   By default, sibling branches in a polytomy that carry identical substitutions are grouped under a new internal node. The step and the flag apply only to the sparse representation.
* `--no-hoist-reversions` — Disable the reversion hoist in polytomies.

   By default, when a child branch carries the exact reversion of a substitution on the node's parent branch, a new node is inserted that groups the node with that child, removing one mutation per reverted position. When set, those reversions are left in place. The step and the flag apply only to the sparse representation.
* `--model-params <MODEL_PARAMS>` — Parameters for the model selected by `--model`, given as space-separated `key=value` items

   The parameters apply to the model of every partition; each `--alignment` file is one partition. In commands with `--partition`, each entry carries its own model and parameters, so `--model-params` is an error next to entries. Keys: `mu` (the overall substitution rate, every named model), `kappa` (the transition to transversion ratio of `k80`, `hky85` and `t92`), `kappa1` and `kappa2` (the rate ratios of `tn93`), `pi` (the equilibrium frequencies of `f81`, `hky85` and `tn93`: comma-separated numbers, one frequency per state of the alphabet in its order, for example A, C, G and T for nucleotides, normalized to sum 1; also `pis`, `Pi`, `Pis`) and `pi_GC` (the GC content of `t92`, between 0 and 1). A key left out keeps its default. Example: `--model hky85 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`.
* `--model <MODEL>` — Substitution model to use

   `infer` infers a model from the data of each partition; each `--alignment` file is one partition. Alternatively, name the model type. If the named model takes parameters, use `--model-params` to set them. The model applies to every file. In commands with `--partition`, each entry carries its own model, so `--model` is an error next to entries.

   [default: infer]

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

* `--partition-rates <PARTITION_RATES>` — Relative substitution rates of the partitions

   Each `--alignment` file is one partition, or each `--partition` entry when there are entries. With two or more partitions, all partitions share the branch lengths, and each partition has one relative rate that scales them for the sites of that partition. `estimate` fits the rate of each partition with the rest of its substitution model fixed (a partition with `--model infer` starts from its fitted rate), then divides all rates by their mean weighted by the number of sites of each partition, so that this mean is 1. The branch lengths are then expected substitutions per site, averaged over all sites of all partitions. `equal` fixes the rate of every partition at 1.

   The setting applies only with two or more partitions.

   [default: estimate]

  Possible values:
  - `estimate`:
    Estimate one relative rate per partition, normalized to a site-weighted mean of 1
  - `equal`:
    Fix the relative rate of every partition at 1

* `--representation <REPRESENTATION>` — Representation of the sequences on the tree

   The dense representation stores the probabilities of every state at every alignment position of each node; the sparse representation stores only the variable positions. Dense is more accurate when branches are long and many sites change, but uses more memory.

   [default: auto]

  Possible values:
  - `auto`:
    Let TreeTime choose the representation; it chooses sparse
  - `dense`:
    Probabilities of every state at every alignment position
  - `sparse`:
    Only the variable positions

* `--gap-fill <GAP_FILL>` — How to handle the gap characters of the input sequences

   `only-terminal` (default, as in TreeTime 0.x): replace the leading and trailing gaps of each sequence with the unknown character of the alphabet (`N` or `X`). `all`: replace every gap with the unknown character. `none`: keep every gap.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `-O`, `--output-all <DIR>` — Folder that receives the default outputs of the command

   Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are created.

   Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs with `{alignment}` in their file name are also written once per file.

   The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}` or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is required when more than one CDS is reconstructed; folder names are used as given.
* `--newick-annotations <NEWICK_ANNOTATIONS>` — Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.

   `plain` writes no comments. `beast` writes a comment such as `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date, and its inferred trait value, as far as the command infers them. Nexus outputs always carry `beast` comments, which FigTree reads.

  Default value: `plain`

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <FILE>` — Path to output Newick tree file.

   `--newick-annotations` sets the node comments.
* `--output-tree-nexus <FILE>` — Path to output Nexus tree file.

   The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
* `--output-tree-auspice <FILE>` — Path to output Auspice v2 JSON tree file
* `--output-tree-mat-pb <FILE>` — Path to output UShER MAT protobuf tree file
* `--output-tree-mat-json <FILE>` — Path to output UShER MAT JSON tree file
* `--output-tree-graph-json <FILE>` — Path to output internal graph JSON tree file
* `--output-tree-dot <FILE>` — Path to output Graphviz DOT tree file
* `-s`, `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`, `node-data`, `model`

* `--no-indels` — Disable indel (insertion and deletion) contributions to branch-length optimization and branch-length distributions.

   When set, branch-length optimization uses the substitution-only likelihood, and the branch distributions of the time inference exclude the Poisson indel term. This matches standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST) and TreeTime 0.x. By default, indels count.
* `--divergence-units <DIVERGENCE_UNITS>` — Units of the divergence values in the augur node data and Auspice outputs.

   `mutations-per-site` (default): branch divergence as substitutions per site. `mutations`: count of reconstructed substitutions per branch that change the state: gaps, `N`, and ambiguity codes compatible with the parent state do not count. `mutations` requires reconstructed sequences, so it does not work with input branch lengths (`--branch-length-mode=input`).

  Default value: `mutations-per-site`

  Possible values: `mutations-per-site`, `mutations`

* `--output-node-data <FILE>` — Path to output augur-compatible node data JSON.

   Contains per-node optimized branch lengths (divergence, substitutions per site) and the input alignment and tree paths. The output is compatible with augur export v2 --node-data.
* `--output-model <FILE>` — Path to output the fitted substitution models as JSON
* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-tree <FILE>` — Tree whose leaf order sets the target order

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-list <FILE>` — File with leaf names in the target order, one name per line

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-aggregate <TOPOLOGY_ORDER_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime prune`

Collapse short branches and branches without mutations, and remove named nodes from a tree

**Usage:**

```sh
treetime prune --tree <FILE> [--alignment <FILE>]... (-O <DIR> | --output-* <FILE>) [OPTIONS]
treetime prune --config <FILE> [OPTIONS]
```

**Examples:**

```sh
treetime prune -t tree.nwk --collapse-short 1e-6 --output-tree-nwk pruned.nwk
treetime prune -t tree.nwk -a aln.fasta.xz --collapse-empty --merge-shared-mutations -O results
treetime prune -t tree.nwk --remove-file outliers.txt --output-tree-nwk pruned.nwk
```

###### **Options:**

* `--config <FILE>` — Config file (YAML or JSON) with the settings of the command

   Command-line flags take precedence over the file. A relative path in the file resolves from the folder of the file, and from the working directory for standard input; a relative path in a flag resolves from the working directory. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `--dry-run` — Print the effective settings as a config file and list the planned output files, then stop without reading inputs or writing files
* `-a`, `--alignment <FILE>` — Aligned FASTA input

   Each `--alignment` file is numbered by its position in the list (1, 2, ...) and is one partition, unless `--partition` entries of `ancestral`, `optimize`, or `timetree` split it into several. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin); `-` may appear at most once.
* `-t`, `--tree <FILE>` — Tree in Newick format

   A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--newick-dialect <NEWICK_DIALECT>` — Newick dialect of the tree inputs (the tree and the reference topology)

   A dialect is a tree structure and an annotation convention, written as `<structure>,<annotations>`. The structure is classic (plain labels and branch lengths), enewick (hybrid nodes such as `x#H1`) or rich (eNewick with support and probability fields). The annotations are plain (comments are text), beast (`[&key=value]`), nhx (`[&&NHX:key=value]`) or mrbayes (`[&B name value]`). For example, `enewick,beast` reads networks with BEAST annotations and `classic,nhx` reads NHX trees. TreeTime reads trees only, so a network input is an error.

  Default value: `classic,beast`

  Possible values: `classic,plain`, `classic,beast`, `classic,nhx`, `classic,mrbayes`, `enewick,plain`, `enewick,beast`, `enewick,nhx`, `enewick,mrbayes`, `rich,plain`, `rich,beast`, `rich,nhx`, `rich,mrbayes`

* `--allow-missing-sequences` — Treat tree tips that have no sequence in the alignment as fully ambiguous (missing data) instead of aborting.

   Without this flag the run aborts when more than one third of the tips lack a sequence, because the alignment may not belong to the tree. Useful for per-CDS translations, where some samples have no peptide for a given CDS.
* `--alphabet <ALPHABET>` — Alphabet of the sequences

   `nuc` (default): nucleotides A, C, G and T, with the IUPAC ambiguity codes and `N` for an unknown base. `aa`: the 20 amino acids and the stop codon `*` as a 21st state, with `X` for an unknown residue. `aa-no-stop`: the 20 amino acids without a state for the stop codon.

  Possible values: `nuc`, `aa`, `aa-no-stop`

* `--zero-based` — Report sequence positions counted from 0 instead of 1.

   Applies to the mutations in the comments of Newick and Nexus trees and to the positions of the homoplasy report. Auspice JSON, augur node data, and UShER MAT outputs count from 1, as their formats require. An input table with positions, such as `--drms`, is always counted from 1.
* `--gap-fill <GAP_FILL>` — How to handle the gap characters of the input sequences

   `only-terminal` (default, as in TreeTime 0.x): replace the leading and trailing gaps of each sequence with the unknown character of the alphabet (`N` or `X`). `all`: replace every gap with the unknown character. `none`: keep every gap.

  Default value: `only-terminal`

  Possible values: `only-terminal`, `all`, `none`

* `--collapse-short <LEN>` — Collapse branches shorter than this length

   Removes each internal branch with a length below LEN and attaches the children of its node to the parent node.
* `--collapse-empty` — Collapse branches without mutations

   Removes each internal branch that has no mutation or other state transition mapped to it.

   Requires --alignment
* `--merge-shared-mutations` — Merge branches in polytomies that share identical mutations

   When sibling branches in a polytomy (node with >2 children) carry identical substitutions, they are grouped under a new internal node. The shared mutations move to the new branch (parent to new node), and only unique mutations remain on children's edges. Reduces tree builder artifacts from arbitrary binary resolution of polytomies.

   Requires --alignment
* `--remove <NAME>` — Names of nodes to remove from the tree

   Accepts several names, separated by commas or spaces.
* `--remove-file <FILE>` — File with names of nodes to remove from the tree, one name per line

   Surrounding whitespace is trimmed and empty lines are ignored. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file. Use `-` to read uncompressed input from standard input (stdin).
* `-O`, `--output-all <DIR>` — Folder that receives the default outputs of the command

   Writes the default outputs of the command as `<dir>/<command>.<ext>`. `--output-selection` restricts the outputs that the folder receives. The per-file flags (`--output-*`) override or supplement the files of the folder.

   An output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails.

   Every output path follows these rules. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` writes a compressed file. The compression level is 2, unless the environment variable `GZ_COMPRESSION`, `BZ2_COMPRESSION`, `XZ_COMPRESSION`, or `ZST_COMPRESSION` sets another level. `.gz`, `.xz`, and `.zst` compress with all `--jobs` threads. `.bz2` compresses with one thread and is much slower to write and to read than the other formats. The path `-` of a per-file flag writes the output uncompressed to standard output. Missing parent folders are created.

   Table outputs choose the delimiter by the extension of their path: `.tsv` writes tabs and `.csv` commas, also before a compression suffix (`dates.tsv.gz`); `-` writes tabs to standard output. Default names use `.tsv`.

   With several `--alignment` files, the outputs that carry sequence positions are written once per file: reconstructed sequences, the augur node data of `ancestral`, Auspice, UShER MAT, Nexus, and Newick with `--newick-annotations` other than `plain`. Files are numbered 1, 2, ... in the order of the `--alignment` list. Their default names carry the number after the command name (`<dir>/<command>.1.<ext>`), and the file name of a per-file flag for such an output must contain `{alignment}`, which is replaced by the number. Other tree outputs with `{alignment}` in their file name are also written once per file.

   The outputs written once per CDS, such as the reconstructed amino-acid sequences, take the placeholder `{cds}` or `%GENE` in the file name of their per-file flag, which is replaced by the name of the CDS. The placeholder is required when more than one CDS is reconstructed; folder names are used as given.
* `--newick-annotations <NEWICK_ANNOTATIONS>` — Node comments of the Newick outputs: `plain`, `beast`, or `nhx`.

   `plain` writes no comments. `beast` writes a comment such as `[&mutations="A45G,G136T",date=2016.45]` after each node, and `nhx` writes the same values as `[&&NHX:key=value]`. The comments hold the mutations on the branch above the node, its date, and its inferred trait value, as far as the command infers them. Nexus outputs always carry `beast` comments, which FigTree reads.

  Default value: `plain`

  Possible values: `plain`, `beast`, `nhx`

* `--output-tree-nwk <FILE>` — Path to output Newick tree file.

   `--newick-annotations` sets the node comments.
* `--output-tree-nexus <FILE>` — Path to output Nexus tree file.

   The nodes carry BEAST comments (`[&mutations="A45G,G136T",date=2016.45]`).
* `--output-tree-auspice <FILE>` — Path to output Auspice v2 JSON tree file
* `--output-tree-mat-pb <FILE>` — Path to output UShER MAT protobuf tree file
* `--output-tree-mat-json <FILE>` — Path to output UShER MAT JSON tree file
* `--output-tree-graph-json <FILE>` — Path to output internal graph JSON tree file
* `--output-tree-dot <FILE>` — Path to output Graphviz DOT tree file
* `-s`, `--output-selection <OUTPUT_SELECTION>` — List of outputs to produce with `--output-all`.

   Restricts which outputs `--output-all` writes. Special value `all` expands to every output available for this command. Requires `--output-all`. Per-file flags are always honored regardless of this selection.

   A selected output that the run has no data for is skipped without a message, for example the substitution model of a run that fits none. A per-file flag for such an output fails. An output that the settings rule out fails in both forms, for example the sequences of `timetree --branch-length-mode input`.

  Possible values: `all`, `nwk`, `nexus`, `auspice`, `mat-pb`, `mat-json`, `graph-json`, `dot`

* `--topology-order <TOPOLOGY_ORDER>` — Canonical topology ordering preset

  Possible values: `keep`, `descendant-count`, `descendant-count-reverse`, `height`, `height-reverse`, `divergence`, `divergence-reverse`, `label`, `label-reverse`, `target-order`, `target-order-reverse`

* `--topology-order-tree <FILE>` — Tree whose leaf order sets the target order

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-list <FILE>` — File with leaf names in the target order, one name per line

   Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`, the leaf order of the input tree is the target order. A path that ends in `.gz`, `.bz2`, `.xz`, or `.zst` reads a compressed file.
* `--topology-order-aggregate <TOPOLOGY_ORDER_AGGREGATE>` — Aggregate used to map a subtree to a target-order position

  Default value: `mean`

  Possible values: `mean`, `median`

* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime pipeline`

Run several analysis commands in order from one config file

The config file (JSON or YAML) lists named steps, each an analysis command with its arguments. The steps run one after another on one dataset, in one process; a step may reference the outputs of an earlier step with `{{ steps.<name>.outputs.<selection> }}`. Use `--dry-run` to print the resolved plan without running anything, and `--steps` to run only some steps.

**Usage:**

```sh
treetime pipeline --config <FILE> [-O <DIR>] [OPTIONS]
```

**Examples:**

```sh
treetime pipeline --config pipeline.yaml -O results
treetime pipeline --config pipeline.yaml --steps clock,timetree -O results
treetime pipeline --config pipeline.yaml --dry-run
```

###### **Options:**

* `--config <FILE>` — Pipeline configuration file (JSON or YAML) describing an ordered list of steps.

   A relative path in the file resolves from the folder of the file, and from the working directory when the file is read from standard input (`-`).
* `-O`, `--output-all <DIR>` — Folder of the outputs: each step writes into a subfolder named after the step. Takes precedence over the top-level `output_all` of the pipeline file. A relative path resolves from the working directory
* `--steps <STEPS>` — Run only these named steps, in the pipeline's list order.

   A referenced upstream step that is not selected must already have its outputs on disk.
* `--dry-run` — Resolve and print the plan (steps, inputs, outputs) without running any step
* `-j`, `--jobs <N>` — Number of threads of the run. Without this flag, the run uses every available CPU thread



## `treetime schema`

Write the JSON schemas of the config files and of the TreeTime data types

**Usage:**

```sh
treetime schema <COMMAND>
```

###### **Subcommands:**

* `write` — Write the JSON schema of a config file, or of every config file and TreeTime data type



## `treetime schema write`

Write the JSON schema of a config file, or of every config file and TreeTime data type

**Usage:**

```sh
treetime schema write [--for <TARGET>] [-o <PATH>] [OPTIONS]
```

**Examples:**

```sh
treetime schema write --for timetree -o timetree.schema.json
treetime schema write --for all -o schemas
```

###### **Options:**

* `--for <TARGET>` — Config file or data type of the schema; `all` writes every schema into the folder of `-o`

  Default value: `all`

  Possible values: `all`, `version-info`, `progress-event`, `error-response`, `pipeline`, `timetree`, `optimize`, `prune`, `ancestral`, `homoplasy`, `clock`, `mugration`

* `-o`, `--output <PATH>` — Output file, or `-` for standard output; with `--for all`, the output folder (default: the working folder)
* `--dry-run` — Print the files that the command would write, then stop without writing them



## `treetime completions`

Print the shell completions of TreeTime

Write the output to the completions file of your shell; the documentation of the shell names its location.

**Usage:**

```sh
treetime completions [SHELL] [OPTIONS]
```

**Examples:**

```sh
treetime completions bash > ~/.local/share/bash-completion/completions/treetime
treetime completions zsh > ~/.zfunc/_treetime
```

###### **Arguments:**

* `<SHELL>` — Name of the shell to generate appropriate completions

  Default value: `bash`

  Possible values: `bash`, `elvish`, `fish`, `fig`, `powershell`, `zsh`




## `treetime help-markdown`

Print the command-line reference in Markdown

**Usage:**

```sh
treetime help-markdown [OPTIONS] > <FILE>
```

**Examples:**

```sh
treetime help-markdown > reference.md
treetime help-markdown | less
```



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

