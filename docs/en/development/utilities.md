---
title: Utilities
description: The shared helpers of the repository, to check before writing a helper, wrapper, serializer, or test scaffold.
---

Check repo utilities before writing helpers, wrappers, serializers, test scaffolding. The index below is updated on a best-effort basis, can get out of date.

- Library functions > project utilities > macros > local one-offs
- Verify a candidate in source (signature, behavior, existing call sites) before use
- Use `treetime-utils` helpers for ndarray, assertions, serialization, and I/O instead of hand-written loops, wrappers, or raw stdlib and third-party calls
- Update this index when adding, moving, or removing a shared utility

## `treetime-utils` index

Crate root: `packages/treetime-utils/src/` (module list in `lib.rs`).

### Errors (`error`)

- Macros: `make_error!()`, `make_report!()`, `make_internal_error!()`, `make_internal_report!()`
- Conversion: `report_to_string()`, `to_eyre_error()`

### Strings and formatting (`fmt`)

- Owned string literal: `o!()`
- `fmt::float`: `float_to_significant_digits()`, `float_to_digits()`, `float_to_exponential()` (C `printf` `%.Ne` form such as `6.617e-02`), `float_to_power_of_ten()` and `float_to_power_of_ten_with_error()` (display form with superscript exponent and Unicode minus, such as `5.07 × 10⁻⁵` and `5.07 ± 0.31 × 10⁻⁵`), extension trait `FloatFormatExt` (`.to_significant_digits()`, `.to_decimal_digits()`, `.to_digits()`)
- `fmt::string`: `quote()`, `quote_single()`, `truncate()` with `TruncateDirection`, `truncate_right_with_ellipsis()`, `name_list()` (first names of a list and the count of the rest, for messages), `counted()` (a count with its noun in singular or plural)

### Collections and iterators (`collections`, `iterator`)

- Vector macros: `vec_of_owned![]`, `vec_u8![]`
- `collections::container`: `get_one()`, `get_exactly_one()`, `minmax()`
- `collections::manyzip`: `Manyzip` (zip any number of iterators into `Vec` items)
- `iterator`: `iterator_intersection()`, `iterator_difference()`, trait `MeanByKey` (`.mean_by_key()`)
- Literal maps and sets: `btreemap!{}`, `btreeset!{}` (workspace dependency `maplit`)
- `IndexSet`, `IndexMap` operations (`.difference()`)

### Ndarray (`array`)

`array::ndarray`:

- Deterministic argmax (Python parity): `argmax_first()`; within a tolerance of the maximum, NaN ignored: `argmax_first_within()` (first index), `argmax_all_within()` (every index)
- Scalar reductions: `max_or()`, `min_or()`, `is_max_above()`
- Elementwise: `minimum()`, `maximum()`, `clamp()`, `clamp_min()`, `clamp_max()`
- Math: `exp()`, `log()`, `cumsum_axis()`
- Creation and padding: `zeros()`
- 1D helpers: `first()`, `last()`, `reverse()`, `reverse_inplace()`, `sorted()`, `sort_inplace()`, `partition_point()` (binary search on an owned array or a view, any layout accepted)
- Matrix construction: `outer()`

Other `array` modules:

- `array::softmax_with_log_norm`: `softmax_with_log_norm_inplace()` (normalizes a mutable view in place and returns the log normalization)

### Serialization adapters (`adapters`)

Use with `#[deser(as = ...)]`:

- `ArrayVec`: `Array1` as a list
- `Array2Rows`: `Array2` as a list of rows
- `TrueOrNull`: a `bool` written as `true` or left out

### Numerics (`least_squares`)

- `LineFit::least_squares()`: ordinary least-squares line fit (`slope`, `intercept`)

### Numeric conversions (`cast`)

Use these instead of `as` casts, which `clippy::as_conversions` rejects:

- `ToF64Exact` (`.to_f64_exact()`): integer to `f64`; lossless through `f64::from` up to 32 bits, checked by a debug assertion to stay within 2^53 for 64-bit types
- `usize_to_u64()`, `u32_to_usize()`: lossless widening, guarded by a compile-time width assertion
- `f64_to_usize_saturating()`: truncation toward zero with the saturation of `as` (NaN and negatives give 0, values above the range give `usize::MAX`)
- `f64_to_i32_saturating()`: the same for `i32` (NaN gives 0, values beyond the range give `i32::MIN` or `i32::MAX`), for pixel coordinates
- `f32` to `f64`: `f64::from`; `u8` to `char`: `char::from`; 128-entry ASCII tables: `AsciiChar::table()` of `treetime-primitives`

### Intervals (`interval`)

- Types: `Range`, `RangeCollection`
- Queries: `range_contains()`, `range_contains_iter()`
- Set operations: `range_union()`, `range_intersection()`, `range_difference()`, `range_complement()`, plus `_iter` variants, `interval_sets_union()`, `interval_sets_intersection()`
- Conversion: `to_interval_set()`, `to_interval_sets()`, `from_interval_set()`

### Dates (`datetime`)

- Parsing: `parse_date()`, `parse_date_with_formats()`, `parse_date_with_format()`, `parse_date_range()`, `parse_date_uncertain()`, `parse_datetime*()`, `DATE_FORMATS`, `DateParserOptions`
- Ranges: `DateRange` (`from_iso()`, `from_ymd()`, `mean()`)
- Year fractions: `date_to_year_fraction()`, `year_fraction_to_date()`, `year_fraction_to_datestring()`, `year_fraction_days_between()`, `date_range_to_year_fraction_range()`, `year_fraction_at_month_start()` (the instant at the start of a calendar month, for date axes)
- Conversion and formatting: `iso()`, `days_in_month()`, `date_from_iso()`, `date_format*()`

### I/O (`io`)

- `io::json`: `json_read_file()`, `json_read_str()`, `json_read()`, `json_write_file()`, `json_write_str()`, `json_write()`, `JsonPretty`; `json_value_read_str()`, `json_value_read_file()`, `to_json_value()`, `from_json_value()` (conversion between deser types and `serde_json::Value`)
- `io::file`: `read_file_with()` and `write_file_with()` (open a path or `-` for standard input or output, decompress or compress by extension, add the path to errors once), `open_file_or_stdin()`, `create_file_or_stdout()`, `FileWriter` (`finish()` returns flush and compression errors), `is_path_stdin()`, `is_path_stdout()`
- `io::fs`: `ensure_dir()`, `absolute_path()`, `extension()`, `read_file_to_string()`, `read_file_to_string_if_exists()` (app state files: `None` when the file is missing)
- `io::compression`: `Decompressor`, `Compressor`, `CompressionType`, `guess_compression_from_filepath()`, `remove_compression_ext()`, `COMPRESSION_EXTENSIONS`
- `io::console`: `is_tty()`
- `env`: `env_var_optional()`

### Concurrency and randomness (`sync`)

- `sync::random`: `get_random_number_generator(seed)` (Xoshiro256++ from a required `u64` seed), `random_choice()`, `random_choice_maybe()`, `random_remove()`

### Initialization (`init`)

- `init::global`: `global_init()` (error report hook of tests and tools), `setup_logger(level, colors)`, `stderr_has_colors()` (standard error is a terminal and `NO_COLOR` is unset), `install_panic_hook(colors)` and `install_error_hook(colors)` (the report hooks of the CLI, without the location and environment sections)
- `init::clap_styles`: `styles()`
- `init::cpu`: `require_cpu_features()` (stop with a message and exit code 1 when the CPU lacks a feature the build enables, called first in the constructors of the binaries), `cpu_has_build_features()`, `CPU_FEATURES_MESSAGE`
- `init::openblas`: `get_openblas_info()`, `get_openblas_info_str()`

### Test assertions (`testing`)

Use these for readable ndarray diffs, float maps with tolerances, proptest failures without backtraces, and exact error messages.

- Scalar and generic: `pretty_assert_eq!()`, `pretty_assert_ulps_eq!()`, `pretty_assert_relative_eq!()`, `pretty_assert_abs_diff_eq!()`, `pretty_assert_neg_inf!()`, `assert_error!()`
- Maps: `pretty_assert_map_ulps_eq!()`, `pretty_assert_map_abs_diff_eq!()`, trait `MapLike` (`BTreeMap`, `HashMap`, `IndexMap`)
- Arrays: `pretty_assert_array_eq!()`, `pretty_assert_array_finite!()`, `pretty_assert_array_nonneg!()`, `pretty_assert_array_positive!()`, `pretty_assert_array_diag_abs!()`, `pretty_assert_array_offdiag_upper_bounded!()`
- Property: `prop_assert_{abs_diff,ulps,relative}_eq!()`, `prop_assert_array_{abs_diff,ulps,relative}_eq!()`, `prop_assert_map_{abs_diff,ulps,relative}_eq!()`
- Property structural: `prop_assert_array_finite!()`, `prop_assert_array_nonneg!()`, `prop_assert_array_positive!()`, `prop_assert_array_upper_bounded!()`, `prop_assert_array_diag_abs!()`, `prop_assert_array_offdiag_nonneg!()`, `prop_assert_array_offdiag_positive!()`, `prop_assert_array_diag_nonpositive!()`
- Tolerances: the ulps and relative assertions compare without an absolute floor unless the call gives `epsilon`, and the absolute-difference assertions require `epsilon`, so a comparison of small values or of a value with 0 states its absolute tolerance at the call site. `clippy.toml` bans the assertions and the comparison macros of `approx` (`ulps_eq!`, `relative_eq!`, `abs_diff_eq!` and their `_ne` forms), whose default `epsilon = f64::EPSILON` accepts any two values closer than about 2.2e-16. A comparison outside an assertion uses the rule of its domain (`coordinates_coincide()`, `YAxisPolicy::same_value()`) or a builder of `approx` (`Ulps`, `Relative`, `AbsDiff`) with an explicit `.epsilon(..)`
- Also check existing support files in `__tests__/` and `pretty_assertions`
- `treetime` crate tests: `RecordingLog` in `packages/treetime/src/test_utils/progress.rs` (a `LogSink` that keeps every message and run warning)

## Other locations

- `packages/treetime-io/src/`: one module per file format: Newick (`nwk`), Nexus (`nex`), FASTA, GFF, Graphviz, Auspice types, UShER MAT, tables (`csv` with `TableFormat` and `CsvWriter`, `dates_csv`, `discrete_states_csv`), name lists (`name_list`), graph writing (`graph`). Entry points follow one naming scheme: `<format>_read(reader)` and `<format>_write(writer)` hold the format logic, `<format>_read_file(path)` and `<format>_write_file(path)` wrap them through `read_file_with()` and `write_file_with()`, and a `_str` form exists only where production code needs a string. `FastaWriter` and `CsvWriter::create()` write records as they arrive and end with `finish()`. PhyloXML is in `packages/util-phyloxml`, not yet wired into `treetime-io`
- `packages/app-output/src/table_output.rs`: `table_write_file()`, `table_create()`, `table_read_file()` (CSV and TSV outputs, with the `TableFormat` of their `OutputSelection`)
- `packages/app-output/src/command_outputs.rs`: `CommandKind` (`outputs()`, the output table of a command with a `DefaultRule` per output; `all_selectable()`, `default_outputs()`, `added_outputs()` for the `OutputSettings` of a run), the one owner of which outputs a command has and writes by default
- `packages/app-output/src/output_conditions.rs`: `OutputSettings` (the settings of a run that decide which outputs exist), `SequenceOutputs`, `NotReconstructed` (`setting()`), `OutputCondition` (`description()`, `requirement()`), the errors `UnavailableOutput` (an output that the run cannot write) and `UnmetCondition` (an output whose condition the run did not meet)
- `packages/app-output/src/output_format.rs`: `table_format()` (TSV or CSV from the file name, TSV for standard output), `FileFormat::check()` (reject an output file name whose extension the output cannot write, with `UnsupportedExtension`), `ChartFormat::from_path()` (SVG, PNG or PDF from the file name, SVG for standard output; the one owner of chart file formats)
- `packages/app-output/src/dates_table.rs`: `write_dates_table_file()` (the node dates table, with the 90% bounds when the run has them)
- `packages/app-output/src/annotated_graph.rs`: `AnnotatedGraph` (the struct of per-node facts that every tree, graph, and per-node table output reads; fact groups `sequences`, `dates`, `traits` are `Option`s), `AnnotatedTreeView::new()` (checked tree handle of it, a `TreeSchedule`, required by the tree-based writers)
- `packages/app-output/src/output_plan.rs`: `ResolvedOutputs::path()` (planned path of a non-tree output), `PlannedFile` with `Requested` (per-file flag or `--output-all`), `output_unavailable()` (error for a named output without data, debug-level skip otherwise), `ResolvedOutputs::ensure_unique_with_expansion()` (destination collision check for outputs expanded into several files), `ResolvedOutputs::alignment_paths()` and `ResolvedOutputs::alignment_tree_outputs()` (the paths of a per-alignment output, by alignment number), `TreeWriteKind::is_per_alignment()` (a tree output written once per alignment file)
- `packages/app-output/src/tree_output.rs`: `write_tree_outputs_per_alignment()` (every per-alignment tree output of a run, one file per alignment, in the coordinates of that file)
- `packages/app-commands/src/commands/shared/output_args.rs`: `OutputRequest` (the destructured output settings of one command, which the conversion of its arguments binds field by field, so a parsed output flag that no code reads fails the lint; `resolve()` builds the output plan)
- `packages/app-commands/src/commands/shared/run_paths.rs`: `RunPaths` (the inputs and outputs of a run by setting and path role, with the translation output template kept apart from the other outputs): `check_without_cds_names()` (reject a run whose output path names one of its inputs, for every path without a CDS template, and a run that reads standard input twice; at argument conversion through `check_command_paths()`), `check_with_cds_names()` (the same comparison after the translation input and output templates are expanded with the CDS names of the run; at argument conversion when `--cds-selection` lists them, and in the `ancestral` run after the CDS selection), `run_inputs()`, `run_outputs()` and `RunPath`; the private `opened_path()` resolves the existing part of a path through symbolic links before the comparison
- `packages/app-cli/src/cli/v0_args.rs`: `rewrite_v0_args()` with the table of TreeTime 0.x spellings (replace or drop a spelling whose meaning is unchanged, stop with the replacement for one whose meaning changed or that v1 removed), the one owner of 0.x spellings
- `packages/app-cli/src/cli/dry_run.rs`: `dry_run_text()` (the `--dry-run` text: the effective config as YAML with resolved paths, then the planned output files with their conditions), using `CommandArgs::planned_paths()` (`packages/app-commands/src/command.rs`, the output paths a run plans to write, per output)
- `packages/treetime-io/src/state_weights_csv.rs`: `state_weights_read_file()` (the state weights table of `mugration --weights`: state name and weight per row)
- `packages/app-commands/src/commands/shared/metadata.rs`: `MetadataArgs` (the `--metadata` flag), `MetadataIdArgs` (`--metadata-id-columns`, `--metadata-delimiters`), shared by every command that reads a metadata table
- `packages/treetime-io/src/csv.rs`: `find_column()` with `ColumnRule::{Exact, OneOf, Containing}` (the column of an input table by an exact requested name, by the first of several candidate names, or by the leftmost header that contains a word; names are normalized like headers), the one owner of column lookup in input tables
- `packages/app-commands/src/commands/shared/sequence_inputs.rs`: `SequenceInputs` (the sequence settings of a command: `read()` reads the alignments, `read_partitions()` reads them and resolves the partitions), `read_tree_and_sources()` (read the tree, every alignment and the partition plans of a sequence command), `read_sequence_sources()` (read each `--alignment` file as one source: apply gap filling, pair the records with leaves, complete missing leaves), `ancestral_input()` (build the `AncestralInput` from the tree and a `SequenceInput`)
- `packages/app-commands/src/commands/shared/partition.rs`: `PartitionArgs` and `PartitionSpec` (the `--partition` flag and config key: one entry per partition with `name`, `alignment`, `sites`, `model`, `model_params`), `partition_entries()` (entries to core `PartitionEntry` values, with the model parameters parsed against each entry's model), `ConfigPartitions` (the same resolution from a config value, for the app run checks), `partition_argument_error()` (a resolver error as the argument error of its settings)
- `packages/treetime/src/partition/resolve.rs`: `resolve_partitions()` (entries, alignment lengths, the default model and the alphabet to `PartitionPlan { layout, model }` values, with the repeated names), `resolve_partition_layouts()` (the site rules alone, giving `PartitionLayout { number, name, alignment, columns }`; prune calls only this step), `check_model_alphabet()` and `model_alphabet()` (a model against the alphabet of a run without alignments), `partition_label()` (`partition <number> '<name>'` in messages), the one owner of partition resolution
- `packages/treetime/src/partition/sites.rs`: `SiteRange` (site ranges `a`, `a-b`, `a-b\k`, 1-based; `FromStr`, `Display`, `columns()` giving 0-based columns; serialized as a string with a schema pattern)
- `packages/treetime/src/partition/split.rs`: `SequenceInput { alignments, partitions }` (the sequence input of the core pipelines: alignment sources and partition plans), `SequenceInput::whole_alignments()` (one plan per alignment with one model, for homoplasy and tests), `split_alignments()` (the node sequences of each partition from its columns; a whole alignment moves without a copy), the one owner of splitting alignments into partitions
- `packages/treetime/src/partition/assembly.rs`: `AlignmentLayout` (whole for an alignment with one partition over all its columns, without per-site data, `is_whole()`; split otherwise, `new()`; `assemble_sequence()`, `assemble_ambiguity()`, `assemble_substitutions()`, `assemble_indels()`, `assemble_edge_mutations()`, `compose_indels()`, `split_indels()`; split pieces carry `InDel.split_from`), `alignment_groups()` (the partitions of a run, given as pairs of an optional `PartitionLayout` and a length, grouped by alignment with their `AlignmentLayout`; build the groups once per topology pass and pass `&[AlignmentGroup]` to each edge operation), the one owner of the mapping between partition and alignment coordinates (`kb/algo/partition-assembly.md`)
- `packages/treetime/src/partition/alignment_sequences.rs`: `stream_alignment_sequences()` (node sequences and mutations of one alignment: a whole alignment streams in one pass, a split alignment is assembled and emitted in preorder), the one owner of streaming per alignment
- `packages/treetime/src/partition/marginal/sparse/alignment_indels.rs`: `sparse_alignment_groups()`, `assembled_edge_indels()`, `composed_edge_pieces()`, `store_edge_indels()`, `merge_sparse_edges()` (assemble, compose and store the indels of sparse edges per alignment), the one owner of indel moves on sparse edges
- `packages/treetime/src/partition/marginal/reconstruction.rs`: `with_layout()` on `DenseReconstruction`, `SparseReconstruction` and `MarginalReconstruction`, `layout()` on `DenseReconstruction`, `SparseReconstruction` and `ReconstructionRef` (the layout of a reconstruction; a reconstruction without one covers its whole alignment); `packages/treetime/src/partition/set.rs`: `PartitionSet::layouts()` and `layouts_or_whole()` (the layout of every partition as a `Cow`, whole layouts filled in), `PartitionSet::alignment_groups()` (the alignment groups of the partitions of a set)
- `packages/treetime/src/partition/result.rs`: `AlignmentResult` (facts of one alignment file: number, path, alphabet, length, mask, indel rate, sequences), `PartitionResult` (facts of one partition: number, name, alignment, sites, model, log-likelihood), `alignment_results()` (the alignment results of a pipeline from its per-alignment outcomes)
- `packages/treetime/src/optimize/partition_rates.rs`: `scale_rates_to_unit_mean()` (divide the rates of partitions by their site-weighted mean; ancestral scales only μ with it, `normalize_partition_rates()` also rescales the branch lengths and indel rates)
- `packages/treetime/src/optimize/gather.rs`: `gather_alignment_edge_indel_counts()` (indel events per edge and alignment, assembled for split alignments, so an indel across partitions counts once)
- `packages/treetime/src/gtr/get_gtr.rs`: `model_states()` (the state characters of a sequence model in the order of `pi` and `W`), `ModelOutput { alignments, models }` (the model document)
- `packages/app-commands/src/commands/shared/ancestral_trees.rs`: `write_ancestral_trees()` (every tree output and the augur node data of a run with reconstructed sequences and branch mutations)
- `packages/treetime/src/homoplasy/classify.rs`: `classify_mutation()` (substitution between canonical states, change involving an ambiguous character, or insertion and deletion), the one owner of this counting rule
- `packages/app-commands/src/commands/shared/gtr_output.rs`: `partition_models()` (the `ModelOutput` of a run from its partition and alignment results), `write_gtr_output()` (the `--output-model` file of a run, or `output_unavailable()` when no model was fitted)
- `packages/treetime-graph/src/key_map.rs`: `NodeMap<V>` and `EdgeMap<V>` (`KeyMap<K, V>` keyed by `GraphNodeKey` or `GraphEdgeKey`), the map of per-node and per-edge values: a vector indexed by key, iterated in ascending key order as `BTreeMap` is, with `entry()`, `retain()`, `keys()`, `Index` that panics on a missing key, `MapLike` for the map assertions, and the `keymap!{}` literal; keep `BTreeMap` only for maps where few nodes have an entry, because the vector holds a slot for every key of the graph
- `packages/treetime-primitives/src/ascii_table.rs`: `AsciiTable<T>` (one value per ASCII character, indexed by `AsciiChar`, sized by `AsciiChar::COUNT`), the one owner of tables indexed by character
- `packages/treetime/src/alphabet/char_slots.rs`: `CharSlots` (`Alphabet::char_slots()`: the slot of each character of an alphabet, for per-character counts such as `Composition`)
- `packages/treetime-graph/src/node_names.rs`: `NodeNames`, the name of every node of a graph, kept next to the graph: `NodeNames::from_input()` (the reader names unnamed input nodes `NODE_<n>` in preorder, renames each later node of a repeated name to `<name>_<k>`, and reports the unnamed nodes and the repeated names), `NodeNames::name_new_node()` (every operation that creates a node names it with the next unused `NODE_<n>`), `names[key]` (the unique name as `&str`), `NodeNames::pairing_name()` (the input name of a renamed node, for matching data from other inputs), the one owner of node names; writers and messages read names and never make them up
- `packages/treetime-graph/src/pair_by_name.rs`: `pair_by_name()` (pair input entries with candidate nodes by their input name: the first entry of a name wins, every node with that input name receives it; returns unmatched entries and repeated entry names), `key_by_name()` (node of each unique name, for names the user types), the one owner of name pairing
- `packages/treetime/src/progress.rs`: `RunWarning`, `RunWarningKind` and `LogSink::warning()` (a typed warning of a run: the CLI and the run log show its message, the run record keeps it for the apps)
- `packages/treetime/src/seq/alignment.rs`: `pair_leaf_sequences()`, `LeafSequences` (alignment records paired with leaves, with the unmatched records for the common length and mask)
- `packages/treetime/src/mugration/trait_pairing.rs`: `pair_trait_rows()` returning `TraitPairing` (metadata trait rows paired with leaves, with the observed values, repeated names, and unmatched rows)
- `packages/treetime/src/ancestral/attach.rs`: `warn_leaves_without_data()` with `LeafData` (warnings for leaves without a sequence or a trait value, which the fit treats as fully ambiguous)
- `packages/treetime/src/gtr/refinement.rs`: `refine_gtr()` with `GtrRefinement` (iteration count, the `InferGtrOptions` of each fit, the optional sampling-bias correction), the one GTR refinement loop, run by mugration on any partition that implements `TransitionCounts` (`packages/treetime/src/partition/marginal/shared/update.rs`, the dense and discrete partitions); `optimize_gtr_rate()` (the rate step, any `MarginalPasses` partition)
- `packages/treetime/src/gtr/infer_gtr/common.rs`: `infer_gtr_from_counts()` with `InferGtrOptions` (GTR fit from transition counts and dwell times), `DEFAULT_PSEUDO_COUNT` (the pseudo-count default of every GTR fit and of `mugration --pc`)
- `packages/treetime/src/alphabet/state_alphabet.rs`: `StateAlphabet::new()` (sorted unique state names of a discrete trait with index lookup; rejects fewer than 2 states)
- `packages/treetime/src/partition/marginal/discrete/state_posterior.rs`: `StatePosterior::new()` (most likely state, profile, and entropy of one node), read by every trait writer
- `packages/treetime/src/partition/marginal/posterior_states.rs`: `most_likely_state()` (the tie rule: first state within 1e-10 of the maximum), `AmbiguityThreshold` with `plausible_states()` (the states of the ambiguity report), the one owner of state decisions over posteriors for traits and sequences
- `packages/treetime/src/seq/sink.rs`: `SeqSink::ambiguity()` (the threshold a sink needs; the stream then computes the ambiguous sites of each emitted node), `SeqItem::seq_with_ambiguity()` (the sequence with the ambiguity codes)
- `packages/app-output/src/ambiguity_outputs.rs`: `AmbiguityOutputs` (`threshold_for()` an output that carries the ambiguity report, `warn_unwritten()`); `packages/app-commands/src/commands/shared/ambiguity.rs`: the `--ambiguity-threshold` and `--ambiguity-outputs` settings of `ancestral`, `timetree`, and `mugration`, resolved with `resolve()`
- `packages/app-commands/src/commands/shared/argument_error.rs`: `argument_error(settings, message)` (a usage error that names the settings at fault; the CLI exits with code 2, and the app marks the settings), `ArgumentError`; `packages/app-commands/src/commands/shared/required.rs`: `missing_required_args()` (the argument error of required flags that neither the command line nor a config file gives, with the flag names from clap), the one owner of argument errors after parsing
- `packages/app-commands/src/commands/shared/input_warnings.rs`: `warn_duplicate_names()` (raise a `RunWarning` for repeated names)
- `packages/app-commands/src/commands/shared/tree_input.rs`: `read_input_tree(path, dialect, log)` (read the input tree of a command in its Newick dialect and warn about repeated node names), `TreeDialectArgs` (the `--newick-dialect <structure>,<annotations>` setting of commands with a tree input, default `classic,beast`), the one owner of the tree input dialect
- `packages/treetime-io/src/nwk.rs`: `graph_from_newick()` (the one conversion from a Newick graph to a TreeTime graph; rejects a network node with more than one parent), `TREE_DIALECT_DEFAULT`
- `packages/app-commands/src/commands/shared/alignment.rs`: `pair_alignment()` (pair FASTA records with leaves at a command input, with FASTA descriptions by node, and warn about repeated sequence names)
- `packages/app-commands/src/commands/shared/dates_input.rs`: `read_input_dates()` (read the dates of a metadata table, pair them with nodes, and warn about repeated rows)
- `packages/treetime-graph/src/schedule.rs`: `GraphSchedule::new()` (validates a rooted graph, also a DAG with several parents and roots, and owns its traversal order; the graph keeps the built index until its topology changes, so build a schedule where a pass needs it instead of passing one through APIs: `preorder()` with every node after all its parents, `postorder()` with every node after all its children, `parents()` and `children()` with edge keys, `roots()`, `leaves()`; errors name the offending node), `require_tree()` (one root and at most one parent per node; algorithms that need a tree call it at their start), `into_tree()` giving `TreeSchedule` (`root()`, `parent()`, `children()`, `preorder()`), the one owner of traversal order. Only the tree-only output writers take `TreeSchedule`; the custom lint `forbidden_item` rejects it in the production code of `treetime`
- `packages/treetime-graph/src/schedule_scan.rs`: `GraphSchedule::scan_down()` and `GraphSchedule::scan_up()` (parallel passes from the roots to the leaves and back: each step gets its node's input and one message per parent or per child, in edge order, and returns its output and its message; no recursion, the first error wins, results do not depend on the thread count), `GraphSchedule::join_parent_edges()` (pair each node's input with one value per parent edge), `ScanOutputs` (`into_map()`, `split_parent_edges()`, `unwrap_shared()` for outputs shared as messages), `PerEdge<T>` (`SmallVec` with inline room for one value, so a tree allocates nothing per node), the one owner of parallel graph passes. A step runs no rayon parallelism of its own, because a waiting `join` inside a step runs stolen work on top of its stack and makes stack use depend on the graph
- `packages/treetime/src/clock/divergence.rs`: `root_to_node_divergences()` (preorder root-to-node sums with a caller-given edge length), `root_to_node_divergences_where_known()` (the same over a `GraphSchedule` that is a tree, `None` below an edge without a length), the one owner of divergence arithmetic
- `packages/treetime-simulation/src/`: generator of test data with a known truth: `simulate_coalescent()` (coalescent tree from `SampleBatch` sampling times and a `TcSchedule`, giving `SimulatedTree` with node dates and edge durations), `TcSchedule` (`constant()`, `piecewise()`, `exponential()`), `SubstitutionModel` (`new()`, `random()`, `jc69()`, `f81()`), `evolve_sequences()` (sequences along the tree, `SimulatedSequences::leaf_records()`, `leaf_records()` for any sequence map), `disjoint_gap_events()` and `apply_gap_events()` (deletions on branches drawn in proportion to their duration, inherited by every node below, never sharing or touching a site), `observe_tips()` with `TipObservation` (leaf sequences as sequenced: IUPAC codes drawn uniformly among the codes with the true base, and runs of N, a coarsening at random), branch lengths `exact()` (expected substitutions per site), `realized()` (fraction of sites that differ between parent and child sequence, gap sites not counted), `perturbed()` (log-normal noise), `FALSE_FAIL_RATE` (the false-failure rate of one simulation test, split across its checks)
- `packages/treetime-graph/src/clades.rs`: `graph_clade_keys()` (the clade of every node as sorted tip ids), `clade_keys()` (the same over nodes in preorder), `first_node_by_key()` (node of each clade), `TipIds` (tip ids shared across graphs), for matching the nodes of two trees by clade
- `packages/treetime/src/test_utils/simulation.rs`: adapters from the generator to the core: `CaseSpec` (`simulate()` gives a `SimulatedCase`), `timetree_cell()` and `timetree_cell_params()` (the standard timetree case), `arb_simulated_case()` (proptest strategy), the input builders `timetree_input()`, `optimize_input()`, `ancestral_input()`, `exact_tip_dates()`, `simulated_sequence_source()`, `first_columns()`, `gtr_from_model()`, `match_internal_nodes()` (true node to inferred node by clade), `run_replicates()` (seeded replicates in parallel), the statistical checks `recovery_check()`, `upper_recovery_check()`, `decrease_check()`, and `require_passed()` (fail with every check's summary)
- `packages/treetime/src/test_utils/reference_log_z.rs`: `reference_log_z()` with `ReferenceInput` (independent oracle of the raw root log mass `log_lh_time` of a fixed tree: per-cell quadrature on uniform grids at h, h/2, h/4 with Richardson extrapolation and an error estimate), `BranchTime::from_likelihoods()`, `NodeFactors` (no prior, or a `CoalescentPrior`); `packages/treetime/src/test_utils/poisson_cherry.rs`: `PoissonCherry` (closed-form log Z, root-age gamma posterior, rate slope of a two-tip tree)
- `packages/treetime/src/test_utils/kingman.rs`: `kingman_neg_log_density()` (Kingman coalescent density of node ages, polytomies as simultaneous mergers), `pair_exposure()`, `node_ages_from_dates()`, `simulated_node_ages()`
- `packages/treetime/src/test_utils/rate_profile.rs`: `FrozenRateState::new()` (branch distributions at the true rate and a point prior at the true node times of a simulated case), `profile()` and `point()` (TreeTime's log Z over clock rates next to `reference_log_z()`)
- `packages/app-backend/src/app_paths.rs`: `app_root()` (app folder of a local app, `TREETIME_APP_DIR`), `AppPaths::resolve()` (profile, runs, logs, and examples folders from the `TREETIME_*_DIR` variables, the `paths:` settings, and the app folder), the one owner of app folder resolution
- `packages/app-backend/src/atomic_write.rs`: `write_atomically()` (write through a temporary file in the same folder, then rename)
- `packages/app-commands/src/config/resolve_paths.rs`: `resolve_config_paths()` (rewrite relative input, input-template, and output paths of a config value to the folder of its config file, using the `x-path` roles of the command schema; `-`, absolute paths, and empty strings stay), `resolve_config_paths_where()` (the same for the settings a predicate selects), the one owner of config path resolution
- `packages/app-commands/src/examples_download.rs`: `examples_url()` (URL of `examples.zip` of this build's release), `download_examples()` (download with progress, unpack into an empty target folder through a private temporary folder); `packages/app-commands/src/examples_unpack.rs`: `unpack_examples()` with `UnpackLimits` (writes only regular files and folders of an archive inside the target, at paths built from checked names; refuses links, special files, duplicate names, mismatched headers, names that leave the folder or that Windows reserves, and archives over the entry or size limit), the one owner of unpacking downloaded archives; `ExampleDownloads` in `packages/app-backend/src/examples_download.rs` (one background download at a time for a local app)
- `packages/app-backend/src/datasets.rs`: `DataIndex::new()` (the index of the example data: every tracked file of `data/` outside hidden files and folders, the datasets with their input files, the example configs and pipelines), written to `data/index.json` by `just data-reindex` (`generate-data-index` of `packages/app-cli`), the one owner of example data discovery; every app reads example data only through it
- `packages/app-ui/build/example-data-vite.ts`: `exampleData()` (the Vite plugin of both web apps that serves `data/index.json` and its files in the dev server and copies them into `dist/data/`). Build code imports generated zod schemas from `@neherlab/app-contracts/zod`, because Vite leaves workspace imports of a config to Node, which cannot load the extensionless entry of the contracts package
- `packages/app-backend/src/runs/lifecycle.rs`: `StartedRun::prepare()`, `compute_run()`, `PreparedRun::complete()` (a run as prepare, compute, and complete steps around the serializable `RunJob` and `RunOutcome`), with `trait RunExecutor` in `packages/app-backend/src/bridge/executor.rs` (where the compute step runs: `ThreadExecutor` in the native apps, a job worker per run in web2), the one owner of the run lifecycle
- `packages/app-ui/src/i18n/locale.tsx` (`@neherlab/app-ui/i18n`): `m` (the compiled Paraglide messages), `LocaleProvider`, `useLocale()`, `locales`, `baseLocale`; every message call passes the locale, `m.<key>(params, { locale })`
- `packages/app-ui/src/chrome/`: `HeaderBar`, `ProductIdentity`, `ThemeCycleButton` (with `nextTheme()`), `SearchTrigger`, the header parts of the app and the docs site
- `packages/design/src/tokens/tokens.ts`: the design tokens (`BASE_COLORS`, `ALIAS_COLORS` with the chart roles, `DATA_COLORS`, `CONTRAST_REQUIREMENTS`, `TYPEFACES`, `CORNERS`, `MOTION`, `resolveColor()`), the one owner of colors, typefaces, and radii; `packages/design/src/tokens/color.ts`: `cssColor()`, `hexColor()`, `hexPair()` (a token as the hex pair of both themes, for Electron windows and `theme-color` tags), `rgbBytes()`, `contrastRatio()`
- `packages/design/src/ui/`: the shared controls (`@neherlab/design/ui/<name>`), with the vendored shadcn/ui components listed in `UPSTREAM.md` and the package's own `wordmark` (product and lab wordmarks), `stat` (a value with its uncertainty), `description-list`, `meter`, `copy-button`, and `data-table`; `packages/design/src/specimen/`: `DesignSpecimen`, the specimen of the docs page
- `packages/app-ui/src/settings/help.tsx`: `SettingSummary` (label slot, flag, default), `SettingHelp` (help and "More" text as Markdown), `defaultText(value, locale)`; `packages/app-ui/src/settings/groups.ts`: `groupedSpecs()` (settings by `--help` group, without the catalog)
- `packages/app-docs/src/lib/urls.ts`: `docsUrl()` and `assetUrl()` (every link of the docs site, with the version base and the locale prefix), `docsPath()` and `assetPath()` for a given base, the one owner of docs URLs
- `packages/app-docs/src/markdown/pipeline.ts`: `docsMarkdownProcessor()` (the Markdown pipeline of the docs site: directives, callouts, steps, math, provided components, version and language URLs)
- `packages/app-docs/src/markdown/provided-imports.ts` with `ProvidedComponents` in `packages/app-docs/src/components/provided-components.ts`: the components that docs pages use without an import, and `Icon<Set><Name>` icons; `remarkProvidedComponents()` adds their imports and `client:` directives, the one owner of components in pages
- `packages/app-docs/src/lib/source.ts`: the page tree of the docs site from the collections `docs` (the `docs/` folder) and `generated`, and the settings pages from the setting catalog
- `packages/app-docs/src/translations/translation-status.ts`: `translationStatus()` and `contentHash()` (missing, stale, and orphan translations of `docs/`, `just docs-i18n-status`)
- `packages/app-ui/src/analysis/inputClaims.ts`: `claimInputs()` (claim input settings before an asynchronous fetch or upload, so a later choice of the user wins); `packages/app-ui/src/analysis/draftInputs.ts`: `applyDraftInputs()` with `assignInputs()` of the draft store (assign inputs and release the uploads they replace), the owners of input choices in the analysis form
- `packages/app-commands/src/json_value.rs`: `JsonValue` (a JSON value without `null`) and `SparseConfig` (settings by key), the named JSON types for free-form fields of schema types; use them instead of `serde_json::Value` in a type that derives `JsonSchema`
- `packages/app-commands/src/command_config.rs`: `CommandConfig` (a command tagged with its parsed config; run records and pipeline steps carry it)
- `packages/app-backend/src/config/choices.rs`: `setting_choices()` (the settings of a command that the form shows as one choice, with the patch each option applies), `active_choices()` (the option a config selects, by the precedence the core uses)
- `packages/treetime-schema/src/no_null.rs`: `NoNull` (schema transform: an `Option` field is optional and never nullable, marked `x-unset`), part of `draft2020_settings()` in `packages/app-commands/src/config/schema.rs`, which every project schema generator uses
- `packages/treetime-distribution/src/distribution_ops/mass_profile.rs`: `MassProfile::new()` and `MassProfile::fractions()` returning `MassFractions` (`mass_between()`, `survival()`, `quantile()`, `survival_inverse()`, `density(x, side)`), the mass model of a time distribution with its tail mass, the one owner of mass, quantile, survival and density queries
- `packages/treetime-grid/src/grid.rs`: `coordinates_coincide()` (two time coordinates are the same point within 4 ulps of a scale), `Grid::coincides()` and `GridFn::coincides()` (the same at the largest magnitude of the grid), the one owner of coordinate comparison
- `packages/treetime-distribution/src/policy.rs`: `YAxisPolicy::same_value()` (two y values are the same probability: relatively for linear probabilities, with an absolute floor at 1 for negative log probabilities)
- `packages/treetime-grid/src/grid_fn.rs`: `GridFn::hard_cell()` (the power-law law of a cell next to a hard end), `GridFn::hard_cell_fit_range()` (the grid points that law is fitted from)
- `packages/treetime-distribution/src/distribution_ops/edge_convolution.rs`: `convolve_forward_across_edge()` and `ForwardConvolution::resolve_posterior()` (forward message on the capped spacing that still resolves the posterior it forms), `ResolvedForward::trim_posterior_to_mass()`
- `packages/treetime/src/timetree/max_posterior_region.rs`: `posterior_survival()` (survival function of a posterior as linear pieces at fixed probability levels `SurvivalLevel`), `max_posterior_region()`
- `packages/treetime-ops/src/convolution.rs`: `convolve_fft()`, `convolve_direct_at()` (one output sample by direct quadrature) and the `ConvolveAlgo` implementations `FftConvolve`, `DirectConvolve`, `NdarrayConvolve` (convolution of two sampled functions, every algorithm by the trapezoid rule over each overlap), the one owner of convolution quadrature
- `packages/treetime-grid/src/max_grid_points.rs`: `MaxGridPoints` (the grid point limit of time inference, validated at least `MaxGridPoints::MIN`), `MaxGridPoints::point_count()` (check a floating-point point count against the limit before any integer conversion or allocation; returns `GridPointLimitExceeded` above it), the one owner of grid size checks
- `packages/app-backend/src/run_limits.rs`: `RunLimits::apply()` (fill an unset `max_grid_points` of a run config from the desktop settings or the server option, and reject values above the server option), the one owner of global run limits
- `packages/app-commands/src/commands/shared/seed.rs`: `SeedArgs` (the `--seed` flag of commands with a random step), `SeedArgs::resolve()` (one seed per run: the given one or a drawn one, logged when the run has a random step), the one owner of run seeds
- `packages/app-commands/src/yaml.rs`: `yaml_read_str()` (typed YAML read with setting paths in errors and YAML 1.1 boolean words), `yaml_value_read_str()` (YAML to `serde_json::Value`, rejects non-finite numbers), `yaml_document()`, `yaml_text()` (YAML output of JSON values through the `saphyr` emitter)
- `packages/treetime/src/lib.rs`: module list of the core library; the workspace `Cargo.toml` lists the crates
- `packages/treetime-openblas`: the OpenBLAS link through pkg-config (`build.rs`, with the install commands of each platform in its error) and the declarations of the `openblas_get_*` functions, which `init::openblas` of `treetime-utils` formats
- `packages/treetime-charts`: `register_chart_fonts()` and `CHART_FONT_FAMILY` (IBM Plex Sans, also registered as `sans-serif`; `FALLBACK_FONT_FAMILY` Lato for the characters that Plex Sans lacks) in `font.rs`, `svg_to_png()`, `svg_to_pdf()` in `convert.rs`, the one owner of chart fonts and SVG conversion; `palette` (generated from the design tokens by `just gen design`: `PAPER`, `INK`, `MUTED`, `GRID`, `AXIS`, `SERIES`, `FIT`, `FAULT`), `ticks` (`linear_ticks()` with one number format per axis, `date_ticks()` at calendar months and years), `text` (`text_width()` as the SVG renderer measures it, `draw_text()` that draws a text in runs by the font that has each character, `Weight`), `scale` (`LinearScale` from data to pixels); compiled only for targets other than WebAssembly. `encode_chart()` in `packages/app-commands/src/rtt_chart/write.rs` writes a chart in the `ChartFormat` of its output file; `RttChartData`, `RttFrame` and `RttSummary` in `packages/app-commands/src/rtt_chart/` hold the points, the axis ranges with the fit and its extrapolation to the root, and the statistics texts that the image chart and the terminal chart share
- `dev/lib/env.sh`: shared shell helpers for scripts in `dev/` (`load_env_file()`, `hash_files()`)
- `dev/lib/dist-flags.sh`: `export_dist_flags <binary|addon> [target]` (RUSTFLAGS and C flags of the shipped builds, with the static GCC unwinder for glibc binaries), the one owner of the shipped CPU flags, used by `dev/cross/build`, the dist, profiling, bench, and addon recipes, `dev/prebuilt-libs-check`, and the build record check of `dev/docker/files/install-image-libraries`
- `kb/algo/` (algorithms, background) and `kb/issues/` (defects, traps): navigation aids, not source substitutes

## Wrong vs correct

```rust
// WRONG - direct constructor noise for static test data
let traits = BTreeMap::from([
  ("A".to_owned(), "usa".to_owned()),
  ("B".to_owned(), "?".to_owned()),
]);

// CORRECT - use project's existing literals
let traits = btreemap! {
  o!("A") => o!("usa"),
  o!("B") => o!("?"),
};
```

```rust
// WRONG - raw serde_json call when repo helper exists
let actual = serde_json::to_string_pretty(&node)?;

// CORRECT - use repository JSON helper
let actual = json_write_str(&node, JsonPretty(true))?;
```

```rust
// WRONG - ad-hoc assertion function
fn assert_profiles_close(lhs: &Array1<f64>, rhs: &Array1<f64>) {
  assert!(lhs.iter().zip(rhs).all(|(a, b)| (a - b).abs() < 1e-10));
}

// CORRECT - shared assertion macro
pretty_assert_abs_diff_eq!(lhs, rhs, epsilon = 1e-10);
```

```rust
// WRONG - custom property assertion wrapper
if !approx::ulps_eq!(lhs, rhs, max_ulps = 4) {
  return Err(TestCaseError::fail("mismatch"));
}

// CORRECT - use shared proptest helper
prop_assert_ulps_eq!(lhs, rhs, max_ulps = 4);
```

```rust
// WRONG - local snapshot wrapper because the adapter was not checked
#[derive(Serialize)]
struct Snapshot {
  profile: Vec<f64>,
}

// CORRECT - put the adapter on the real struct when that is the real contract
#[derive(Serialize, Deserialize)]
struct NodeData {
  #[deser(as = ArrayVec)]
  profile: Array1<f64>,
}
```
