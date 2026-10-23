# Partitions: sets of sites of the alignment files, one model and one relative rate per partition

Ancestral, optimize, timetree and prune accept several `--alignment` files. Ancestral, optimize and timetree also accept `--partition` entries, which split the sites of the files into partitions. A partition is a set of sites (alignment columns) of one file, analyzed with its own substitution model. Without entries, each file is one partition. All partitions of a run share one tree. In optimize and timetree they also share the branch lengths, and each partition has one relative rate that scales the branch lengths for its sites. Indel events and indel rates belong to the file, and every output that carries sequence positions stays one file per alignment, in the coordinates of that file.

Homoplasy stops with an error when it gets more than one file, and clock accepts `--alignment` but does not read the sequences. Prune, homoplasy and clock have no `--partition` setting.

**Type**: Scientific and large design choice for v1.

**Status**: Approved on 2026-10-18 (one partition per file); partitions as sets of sites approved on 2026-10-22.

## Rule

### Alignments and partitions

- **Alignment**: one `--alignment` file. TreeTime numbers the alignments 1, 2, ... in the order of the `--alignment` list. The number identifies the file; the path is a label for outputs and messages, and no code identifies a file by its path or name. Segment builds often use one file name in several folders, for example `results/ha/aligned.fasta` and `results/na/aligned.fasta`, so a file name cannot identify a file. [kb/decisions/duplicate-names-warned-ids-from-input-order.md](duplicate-names-warned-ids-from-input-order.md) states the same principle for all input
- **Partition**: a set of sites of one alignment with its own substitution model: exchange rates W, equilibrium frequencies π and overall rate μ. TreeTime numbers the partitions 1, 2, ... in the order of the `--partition` entries, or one per alignment in alignment order when there are no entries. A name is optional; messages write `partition <number>`, followed by the name in quotes when there is one. Repeated names raise one warning, because the number identifies the partition
- **One representation per run**: all partitions of a run are dense or all are sparse (`enum PartitionSet` in [packages/treetime/src/partition/set.rs#L70](../../packages/treetime/src/partition/set.rs#L70)). `--representation` selects the representation for the whole run
- **Commands that need one file**: homoplasy, `--translations` and `--annotation` give positions in one alignment, so they stop with an error when the run has several files. Several partitions of that one file are allowed, because the outputs stay in the coordinates of the file

### Partition entries

An entry of `--partition` (config key `partition`, one flag per entry) is a YAML or JSON map with the keys `name`, `alignment`, `sites`, `model` and `model_params` (`struct PartitionSpec`, [packages/app-commands/src/commands/shared/partition.rs#L141](../../packages/app-commands/src/commands/shared/partition.rs#L141)). For example, a file with a coding region from site 101 to 1600:

```yaml
alignment:
  - genome.fasta.xz
partition:
  - {name: cp1, sites: ['101-1600\3']}
  - {name: cp2, sites: ['102-1600\3']}
  - {name: cp3, sites: ['103-1600\3'], model: hky85, model_params: ['kappa=4']}
  - {name: noncoding}
```

**Site ranges** (`struct SiteRange`, [packages/treetime/src/partition/sites.rs#L14](../../packages/treetime/src/partition/sites.rs#L14)) are numbered from 1: `a` (one site), `a-b` (the sites a to b, both included) and `a-b\k` (every k-th site from a to b, starting at a). This is the syntax of RAxML-NG partition files and NEXUS charsets, which IQ-TREE, MrBayes and PartitionFinder also read ([raxml-ng src/PartitionInfo.cpp#L48-L79](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/PartitionInfo.cpp#L48-L79)). A range with position 0, a first site after the last site, or a step of 0 is an error. YAML reads a backslash in a double-quoted string as an escape, so ranges are written in single quotes or without quotes.

**Resolution** (`fn resolve_partition_layouts()` and `fn resolve_partitions()`, [packages/treetime/src/partition/resolve.rs#L10-L100](../../packages/treetime/src/partition/resolve.rs#L10-L100)) turns the entries and the alignment lengths into the partitions of the run. The CLI and the apps call the same functions, so they apply the same rules and messages:

- **No entries**: one partition per alignment with all its sites and the model of `--model` and `--model-params`. An alignment without sites is an error
- **Entries without alignment files**: an error, also when `--alignment` is absent from both the flags and the config file
- **Alignment of an entry**: `alignment` must be a number of the `--alignment` list. It is required when the list has two or more files, and every alignment needs at least one entry
- **Sites inside the alignment**: every site range ends at or before the last site of its alignment
- **No site twice**: a site listed twice in one entry, or in two entries of one alignment, is an error that names the partitions and the first repeated site, as RAxML-NG stops on a site that another range already took ([raxml-ng src/PartitionInfo.cpp#L69-L79](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/PartitionInfo.cpp#L69-L79))
- **Remaining sites**: an entry without `sites` takes every site of its alignment that the other entries leave. At most one entry per alignment can leave out `sites`, and such an entry must get at least one site
- **Every site covered**: sites that no entry takes are an error that gives their count and the first of them, because a dropped site changes the total length and the branch-length floor. RAxML-NG also stops on unassigned sites unless its safety check `msa_extra_cols` is off ([raxml-ng src/main.cpp#L469-L497](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/main.cpp#L469-L497))
- **No empty partition**: an entry with an empty `sites` list is an error

The errors reach the user as argument errors (exit code 2) that name the settings at fault: `partition` and `alignment` for the site rules, `partition` and `alphabet` (or `model`, `model_params` and `alphabet` without entries) for the model check. `--dry-run` prints the entries as written, because their resolution needs the alignment lengths. Prune resolves only the layouts, with no entries and no model, and keeps its internal JC69 for its Fitch reconstructions, so an amino-acid prune run stays valid.

### Models

- **One owner of each partition's model**: without entries, `--model` and `--model-params` set the model of every partition. With entries, each entry carries its own `model` (default `infer`, the default of `--model`) and `model_params` (the `key=value` items of `--model-params`, parsed against the entry's model). `--partition` conflicts with `--model` and `--model-params` through one clap rule (`conflicts_with_all`, [packages/app-commands/src/commands/shared/partition.rs#L127](../../packages/app-commands/src/commands/shared/partition.rs#L127)), so the command line, the config check and the conflicts of the app's setting catalog read the same rule. Parameters on an entry whose model is `infer` are an argument error that names the entry
- **Model and alphabet**: every partition's model must fit the alphabet of the run (`--alphabet`, default `nuc`). The models with transition and transversion rates (`k80`, `hky85`, `t92`, `tn93`) need `nuc`, because they set rates between nucleotides. `jtt92` needs `aa-no-stop`, because its matrix holds the 20 amino acids without the stop `*`, and `aa` has 21 states. `jc69` and `f81` build one state per canonical character of the alphabet, as 0.x does ([packages/legacy/treetime/treetime/nuc_models.py](../../packages/legacy/treetime/treetime/nuc_models.py)), so they fit every alphabet; `infer` fits any. Frequencies given with `pi` must have one value per state of the alphabet. The check also covers runs without entries, so `--alphabet aa --model hky85` and `--alphabet aa --model jtt92` stop with an error, while `--alphabet aa --model jc69` runs. A timetree run without `--alignment` checks the model against `--alphabet`, or against the alphabet that the model implies when `--alphabet` is absent (`fn check_model_alphabet()`, [packages/treetime/src/partition/resolve.rs#L44-L50](../../packages/treetime/src/partition/resolve.rs#L44-L50))
- **Entries without effect**: ancestral with `--method-anc parsimony` uses no model and treats every site alone, so split and unsplit runs give the same result; `--partition` there is an argument error that names `partition` and `method_anc`. In timetree, the entries keep their effect when `--branch-length-mode` resolves to `input`: input mode builds no partition, but the models of the partitions set the saturation level of its branch-length likelihoods (see Relative rates)

### Total length and branch-length floor

The total length of a run is the number of sites of all files together, and "one mutation" is 1 / total length. A partition floors a short branch at a fraction of one mutation, so every partition of a file uses the length of the whole file or run, never its own number of sites; otherwise a file split into three partitions would get a floor three times higher:

- **Optimize and timetree**: every partition uses the total length of the run (`fn PartitionSet::total_length()`, [packages/treetime/src/partition/set.rs#L107-L109](../../packages/treetime/src/partition/set.rs#L107-L109)), so every partition clamps a shared branch to the same value
- **Ancestral**: every partition uses the length of its alignment ([packages/treetime/src/ancestral/plan.rs#L124-L168](../../packages/treetime/src/ancestral/plan.rs#L124-L168)), because the partitions of different files share no estimated branch lengths. A split file keeps the floor of the unsplit file

### Relative rates

The switch `--partition-rates` (`struct PartitionRatesArgs`, [packages/app-commands/src/commands/shared/model.rs#L220](../../packages/app-commands/src/commands/shared/model.rs#L220)) applies to optimize and timetree with two or more partitions. The partitions are counted as the entries when there are entries, otherwise as the `--alignment` files, so one file split into codon positions counts as several partitions; with fewer than two, the setting is an argument error (`fn check_partition_rates()`, [packages/app-commands/src/commands/shared/model.rs#L243-L270](../../packages/app-commands/src/commands/shared/model.rs#L243-L270)).

- **`estimate`** (default): the run estimates μ_p of each partition with W and π of the partition fixed, starting from the current μ_p (the fitted μ for `infer`, the μ of the model otherwise). Partition p then expects μ_p * b substitutions per site on a branch of length b
- **`equal`**: every μ_p is 1. An `infer` partition sets μ = 1 after the fit. This setting suits tests and partitions that are known to evolve at the same speed

**Normalization**: after each rate step, `fn normalize_partition_rates()` ([packages/treetime/src/optimize/partition_rates.rs#L70-L88](../../packages/treetime/src/optimize/partition_rates.rs#L70-L88)) computes the site-weighted mean m of the rates, divides every μ_p and every indel rate by m, and multiplies the branch lengths by m. A branch length is then the expected number of substitutions per site, averaged over all sites of all files, and the clock rate has the same unit. The clock rate of partition p is the clock rate of the run times μ_p; the model document reports it per partition.

The normalization does not change the likelihood on branches at or above the floor. A branch below the floor, for example a zero-length branch, enters the propagation as μ_p * floor before the normalization and as (μ_p / m) * floor after it, so its term changes. The floor stays on the branch length, because single-partition runs and v0 floor the branch length. The optimize loop therefore computes the likelihood again after each normalization. When the rates converge, m approaches 1 and the change vanishes.

**Rate step**: `fn optimize_gtr_rate()` ([packages/treetime/src/gtr/refinement.rs#L126](../../packages/treetime/src/gtr/refinement.rs#L126)) evaluates the likelihood at 0.01 * sqrt(μ), sqrt(μ) and 100 * sqrt(μ), and minimizes with a bracketed Brent search in sqrt(μ) when the middle point has the lowest cost. When it does not, the partition-rate step takes the end of the range with the lower cost (`RateRangeEdge::LowerCostEnd`), so a partition whose likelihood keeps rising toward one end moves to that end. After the step, the caller runs the marginal update with the new rates, so the messages always belong to the current rates.

**Schedule**:

- **Optimize**: each iteration runs the branch-length step, the damping, the rate step and the normalization (`fn update_partition_rates()`, [packages/treetime/src/optimize/partition_rates.rs#L14-L32](../../packages/treetime/src/optimize/partition_rates.rs#L14-L32), called from [packages/treetime/src/optimize/run_loop.rs#L76](../../packages/treetime/src/optimize/run_loop.rs#L76)). The stop rule compares likelihoods after the normalization and needs no separate rate tolerance. The best iteration that the loop returns holds its branch lengths, rates and indel rates together
- **Timetree**: each of the two ML branch-length passes before the dating rounds ends with a rate step and a normalization ([packages/treetime/src/timetree/pre_loop.rs#L94](../../packages/treetime/src/timetree/pre_loop.rs#L94)). The dating rounds keep the rates fixed
- **Ancestral** keeps the input branch lengths, and **prune** counts parsimony substitutions only, so neither estimates rates

**Rates after the first fit** (`enum InferredRate`, [packages/treetime/src/partition/inferred_rate.rs](../../packages/treetime/src/partition/inferred_rate.rs)):

- **One partition**: a run with one partition has no rate step and no normalization. An `infer` model gets μ = 1 after the fit, as in v0 ([packages/legacy/treetime/treetime/treeanc.py#L1626-L1631](../../packages/legacy/treetime/treetime/treeanc.py#L1626-L1631)), in dense and sparse mode, so the branch lengths and the clock rate are expected substitutions per site. The rule keeps the branch lengths, as v0 does: optimize and timetree start the branch-length optimization from the input branch lengths. A named model keeps its μ, which is 1 unless `mu` sets it, so a user-given `mu` keeps its meaning
- **Optimize and timetree with several partitions**: each `infer` partition keeps its fitted μ as the start of the rate step. Optimize also normalizes all rates once before its first iteration when any partition has the model `infer` ([packages/treetime/src/optimize/pipeline.rs#L168](../../packages/treetime/src/optimize/pipeline.rs#L168) and [#L186-L188](../../packages/treetime/src/optimize/pipeline.rs#L186-L188)), the rates of named models included; without an `infer` partition, the named models keep their μ
- **Ancestral**: each file is reconstructed on the input branch lengths. A file with one partition applies the one-partition rule, so a run with several files gives the models of separate runs. In a file with two or more partitions, the `infer` partitions keep their fitted μ, scaled together to a site-weighted mean of 1 over those partitions (`fn scale_rates_to_unit_mean()`, [packages/treetime/src/optimize/partition_rates.rs#L90](../../packages/treetime/src/optimize/partition_rates.rs#L90)); partitions with a named model keep the μ of their model. With named models at their default μ of 1, the site-weighted mean of the file is then 1, so the input branch lengths keep their unit, and the third codon position can be faster than the first. This scaling changes only μ, unlike the normalization of optimize and timetree, which also rescales the branch lengths
- **Mugration** keeps its fitted and optimized rate, as v0 does ([packages/legacy/treetime/treetime/wrappers.py#L787-L808](../../packages/legacy/treetime/treetime/wrappers.py#L787-L808))

**Timetree input mode** (`--branch-length-mode input`, or `auto` resolved to it) builds no partition, but takes the saturation level `p0 = 1 - Σ π_i²` of its branch-length likelihoods from the models of the run: the site-weighted mean over the partitions, each partition's `p0` from its own model fitted on its sites (`fn partitions_saturated_difference()`, [packages/treetime/src/timetree/inference/input_branch_length.rs#L262-L288](../../packages/treetime/src/timetree/inference/input_branch_length.rs#L262-L288)). Without entries, this is the mean over the files. A run without `--alignment` takes `p0` from the model of `--model`, `--model-params` and `--alphabet`.

### Indels

- **Events per file**: indel events belong to the alignment file. For a file with one partition, the events are those of the partition. For a split file, the partition pieces of an edge are assembled into the events of the file ([kb/algo/partition-assembly.md](../algo/partition-assembly.md)), so a deletion across the sites of several partitions counts once, and a three-site deletion of a file split into codon positions is one event, not three (`fn gather_alignment_edge_indel_counts()`, [packages/treetime/src/optimize/gather.rs#L98-L120](../../packages/treetime/src/optimize/gather.rs#L98-L120))
- **One indel rate per file**: the indel rate λ_a of file a is the number of indel events of the file divided by the total branch length (`fn IndelPolicy::edge_indels()`, [packages/treetime/src/optimize/indel.rs#L23-L45](../../packages/treetime/src/optimize/indel.rs#L23-L45)). The model document reports it per file
- **One combined term**: all files share the branch lengths, so the sum of the per-file Poisson terms of a branch depends on the branch length as one Poisson term with the summed count and the summed rate Σ_a λ_a. TreeTime evaluates that one term, so every consumer of the indel term of [kb/decisions/optimize-indel-contribution-to-likelihood.md](optimize-indel-contribution-to-likelihood.md) works unchanged
- **File boundaries**: a deletion that crosses the boundary between two files counts as one event in each file. Terminal gap filling at the edge of a file can turn such gaps into unknown characters, so two files cut from one alignment give the results of the unsplit alignment only when no gap touches the cut. Partitions inside one file have no such boundary, because gap filling acts on the whole file before it is split
- **Sparse topology operations**: sibling merges, the reversion hoist, collapses, edge merges and reroots decide on the assembled indels of each file, in file coordinates, and split the result back to the partitions ([kb/algo/partition-assembly.md](../algo/partition-assembly.md)). Partition coordinates would give other decisions: a three-site deletion is one event unsplit and three in codon partitions, and sibling deletions at (0,4) and (0,3) differ unsplit but are equal in the second codon partition. For substitutions, partition and file coordinates give the same decisions, so substitutions keep their per-partition keys

### Outputs

- **Per-file outputs**: reconstructed sequences, augur node data, Auspice JSON, annotated trees and MAT files are one file per alignment, in the coordinates of that file: `<command>.<n>.<ext>` under `--output-all`, or the `{alignment}` placeholder in a per-file flag. A split file streams its assembled node sequences, so sequences, mutations and ambiguity codes equal those of the unsplit file on the same tree
- **Model document** (`--output-model`, `struct ModelOutput`, [packages/treetime/src/gtr/get_gtr.rs#L37-L46](../../packages/treetime/src/gtr/get_gtr.rs#L37-L46)): `alignments` lists each file with its number, path, number of sites and indel rate; `models` holds one entry per partition with its number, name, alignment number, number of sites, clock rate (timetree), state characters and model parameters

## Why

- **Segmented genomes and genes**: influenza segments, or genes with different composition, evolve under different substitution models but share one tree. One file per segment is the input that segment builds already produce
- **Sites inside a file**: codon positions, genes or genome regions inside one alignment evolve under different models and at different speeds. Entries over sites let users fit them without cutting the file, with the syntax of RAxML-NG partition files and NEXUS charsets
- **Every site covered**: a site that no partition takes would silently leave the analysis and change the total length, the branch-length floor and the unit of the branch lengths
- **One owner of the model**: with entries, `--model` and an entry's `model` would both describe the model of one partition, so they cannot be given together
- **Model and alphabet**: `k80`, `hky85`, `t92` and `tn93` name nucleotide pairs, so with an amino-acid alphabet `t92` and `tn93` fail and `k80` and `hky85` build a meaningless model; `jtt92` has 20 states and cannot serve the 21 states of `aa`
- **Relative rates**: partitions evolve at different speeds. With shared branch lengths and one fixed μ per partition, a slow partition would pull the branch lengths down and a fast partition would push them up. One relative rate per partition with shared branch lengths is the default of RAxML-NG (`--brlen scaled`) and the `-p` option of IQ-TREE
- **Normalization**: the likelihood depends on μ_p and b only through the products μ_p * b, so the rates and the branch lengths have one free common scale. The site-weighted mean of 1 fixes that scale and gives branch lengths in substitutions per site averaged over all sites, the unit of a single-partition run with μ = 1
- **Ancestral rates**: ancestral keeps the input branch lengths, so it cannot rescale them. Scaling the inferred rates of a split file to a mean of 1 keeps the unit of the input branch lengths and still lets the partitions of the file differ
- **Range end**: a partition without any change prefers the lowest rate, so its optimum lies outside every range. Keeping the start rate in that case leaves a rate that the data contradict. When the middle point is not the lowest, the end with the lower cost costs at most as much as the start, so the step never lowers the likelihood
- **Indels per file**: an indel removes or adds consecutive sites of the sequence, independent of how the sites are grouped into models. Counting it once per file keeps the event counts, the indel rates and the topology decisions of a split file equal to those of the unsplit file, and one combined term keeps one likelihood formula for all branch-length consumers

## Verification

v0 accepts one alignment and has no counterpart, so the evidence comes from equalities, invariants and exact errors:

- **Split into codon positions**: one alignment as one partition, and the same alignment as three interleaved partitions (`1-L\3`, `2-L\3`, `3-L\3`), with `jc69` and `--partition-rates equal`, give the same results at `1e-8` on Zika20, Flu20 and Ebola20 (the gapped dataset), in dense and sparse mode, with the indel term on:
  - optimize, with the sparse topology operations on: the same branch lengths, edge mutations, sequences and indel rates (`test_optimize_codon_positions_give_the_branch_lengths_and_mutations_of_one_alignment` in `packages/treetime/src/optimize/__tests__/test_pipeline_partitions.rs`)
  - timetree: the same dates, mutations and indel rates (`test_timetree_codon_positions_give_the_dates_and_mutations_of_one_alignment` in `packages/treetime/src/timetree/__tests__/test_pipeline_partitions.rs`)
  - ancestral: the same sequences, mutations, ambiguity codes and total log-likelihood (`test_ancestral_codon_positions_give_the_sequences_and_mutations_of_one_alignment` in `packages/treetime/src/ancestral/__tests__/test_pipeline_split_alignment.rs`)
  - dense cases use the leading 399 sites (optimize) or 198 sites (timetree) for speed; sparse cases use every whole codon of the alignment
- **Models per partition**: an ancestral run with `jc69` on one codon position and `hky85` with κ = 4 on the other two equals two separate runs of the two column subsets with those models, with the separate runs floored at the full alignment length (`test_ancestral_split_alignment_with_two_models_gives_each_partition_the_result_of_its_separate_run`, one case per representation and partition, and `test_ancestral_split_alignment_with_two_models_sums_the_log_likelihoods_of_the_separate_runs`, one case per representation)
- **Floor and ancestral rates**: every partition of a split ancestral file floors a zero-length branch at the alignment length; the inferred rates of a split file have a site-weighted mean of 1, and named partitions keep the rate of their model (`packages/treetime/src/ancestral/__tests__/test_partition_plans.rs`)
- **Optimize rates**: one `infer` and one named partition are normalized before the first iteration; without an `infer` partition the named models keep their rates; one partition keeps the `mu` of its named model (`packages/treetime/src/optimize/__tests__/test_pipeline_partition_models.rs`)
- **Timetree input mode**: the saturation level is the site-weighted mean over the partition models (`test_input_branch_model_takes_the_site_weighted_saturation_of_the_partition_models` in `packages/treetime/src/timetree/inference/__tests__/test_input_branch_length.rs`)
- **Indels**: a deletion across the cut of two partitions counts once per alignment with the summed rate (`test_optimize_deletion_across_the_cut_counts_once_per_alignment_with_the_summed_rate`); each sparse topology operation on a codon-split alignment gives the result of the unsplit alignment, compared by edge name: branch lengths, substitutions, indels in their order, indel event counts, the assembled root sequence, and the node sequences after a marginal update (`packages/treetime/src/optimize/topology/__tests__/test_split_alignment_topology.rs`):
  - sibling merge with a shared substitution and a private three-site deletion (`test_sibling_merge_of_a_codon_split_alignment_gives_the_branch_lengths_of_the_unsplit_alignment`)
  - collapse, which keeps adjacent insertions of two edges as two events (`test_collapse_in_a_codon_split_alignment_keeps_adjacent_insertions_of_two_edges_as_two_events`)
  - reroot with an edge split, edge inversion and the merge of the trivial old root (`test_reroot_of_a_codon_split_alignment_inverts_and_merges_the_events_of_the_unsplit_alignment`)
  - merge of the trivial root with adjacent events on the two merged edges (`test_trivial_root_merge_of_a_codon_split_alignment_composes_the_events_of_the_unsplit_alignment`)
  - the edge merge of `PartitionSet::reconcile_topology()` (`test_edge_merge_of_a_codon_split_alignment_composes_the_events_of_the_unsplit_alignment`)
  - reversion hoist, which tests interaction on alignment positions: deletions adjacent only inside one partition, and adjacent insertions of two partitions (`test_reversion_hoist_of_a_codon_split_alignment_tests_interaction_on_alignment_positions`)
- **Assembly**: the assembly invariants are unit and property tests in `packages/treetime/src/partition/__tests__/test_assembly.rs`
- **Resolution and site ranges**: one test per rule with its exact message, the model and alphabet errors, and a property test that `SiteRange` round-trips through its text (`packages/treetime/src/partition/__tests__/test_resolve.rs`)
- **Several files**: one alignment and the same alignment cut into two column blocks (two files) give the same branch lengths, dates and pruned trees with `jc69` and `--partition-rates equal` when no gap touches the cut. The most likely results do not depend on the order of the files
- **Ancestral factorization**: for the most likely reconstruction, an ancestral run with two files equals two single-file runs, and its log-likelihood equals their sum. Posterior sampling (`--sample-from-profile` with `root` or `all`) draws from one random generator per run, partition after partition in partition order. A sampled result therefore depends on the order of the partitions, and every partition after the first draws other values than in a separate run with the same seed
- **Scale**: the normalized rates and branch lengths do not depend on the scale of the starting rates, and one normalization keeps the likelihood on branches above the floor
- **Rate recovery**: each replicate draws three files of 3000 sites with `treetime-simulation` on the tree and branch lengths of `data/zika/86/tree.nwk`, under JC69 with the rates 0.3, 1.0 and 2.0. The error of a file is the log ratio of its estimated rate and its normalized rate (3/11, 10/11, 20/11). The test passes when the `1 - 2 alpha` confidence interval of the mean error lies within ±0.03 for every file, with `alpha = 1e-3 / 3` per file (two one-sided t-tests). A pilot of 200 replicates (mean errors -0.005, -0.009 and 0.002, standard deviations 0.130, 0.060 and 0.034) sets the size to 10275 replicates, so an estimator within the upper bounds of the pilot fails with a probability of at most about `1e-3`. The test runs in the slow tier, `just test-slow-rs` (`test_sim_optimize_partition_rates_recover_the_simulated_relative_rates` in `packages/treetime/src/optimize/__tests__/test_pipeline_partition_rate_recovery.rs`)

## v0 behavior

v0 reads one alignment into `TreeAnc` with one GTR model ([packages/legacy/treetime/treetime/treeanc.py#L171-L184](../../packages/legacy/treetime/treetime/treeanc.py#L171-L184)). Users concatenate segments into one alignment, which gives every segment and every codon position the same model and the same rate.

## Related

- [kb/decisions/partition-system-architecture.md](partition-system-architecture.md)
- [kb/algo/partition-assembly.md](../algo/partition-assembly.md)
- [kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md)
- [kb/decisions/mugration-rate-step-without-bracket-keeps-start-messages.md](mugration-rate-step-without-bracket-keeps-start-messages.md)
- [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md)
