# Several alignment files: one partition per file, one relative rate per file

Ancestral, optimize, timetree and prune accept several `--alignment` files. Homoplasy stops with an error when it gets more than one file, and clock accepts `--alignment` but does not read the sequences. Each file is one partition with its own substitution model, and all partitions of a run share one tree. In optimize and timetree they also share the branch lengths, and each file has one relative rate that scales the branch lengths for its sites.

**Type**: Scientific and large design choice for v1.

**Status**: Approved on 2026-10-18.

## Rule

### Partitions

- **One partition per file**: each `--alignment` file is one partition and covers all sites of that file. A partition has its own substitution model: exchange rates W, equilibrium frequencies π and overall rate μ. The `--model` choice applies to every file, and `--model infer` fits each file on its own
- **Identity by position**: TreeTime numbers the partitions 1, 2, ... in the order of the `--alignment` list. The file path is metadata: outputs and the apps show it, and no code identifies a file by its path or name. Segment builds often use one file name in several folders, for example `results/ha/aligned.fasta` and `results/na/aligned.fasta`, so a file name cannot identify a file. [kb/decisions/duplicate-names-warned-ids-from-input-order.md](duplicate-names-warned-ids-from-input-order.md) states the same principle for all input
- **One representation per run**: all partitions of a run are dense or all are sparse (`enum PartitionSet` in [packages/treetime/src/partition/set.rs#L49](../../packages/treetime/src/partition/set.rs#L49)). `--dense` selects the representation for the whole run
- **Total length**: the total length of a run is the number of sites of all files together, and "one mutation" is 1 / total length. In optimize and timetree, the floors of short branches use the total length (`fn PartitionSet::total_length()`, [packages/treetime/src/partition/set.rs#L86-L88](../../packages/treetime/src/partition/set.rs#L86-L88)), so every partition clamps a shared branch to the same value. Ancestral keeps the floor of each file, because its partitions share no estimated branch lengths
- **Commands that need one file**: homoplasy, `--translations` and `--annotation` give positions in one alignment, so they stop with an error when the run has several files

### Relative rates

The switch `--partition-rates` (`struct PartitionRatesArgs`, [packages/app-commands/src/commands/shared/model.rs#L84-L101](../../packages/app-commands/src/commands/shared/model.rs#L84-L101)) applies to optimize and timetree with two or more files:

- **`estimate`** (default): the run estimates μ_p of each file with W and π of the file fixed, starting from the current μ_p (the fitted μ for `--model infer`). Partition p then expects μ_p * b substitutions per site on a branch of length b
- **`equal`**: every μ_p is 1. A `--model infer` file sets μ = 1 after the fit. This setting suits tests and files that are known to evolve at the same speed

**Normalization**: after each rate step, `fn normalize_partition_rates()` ([packages/treetime/src/optimize/partition_rates.rs#L69-L104](../../packages/treetime/src/optimize/partition_rates.rs#L69-L104)) computes the site-weighted mean m of the rates, divides every μ_p and every indel rate by m, and multiplies the branch lengths by m. A branch length is then the expected number of substitutions per site, averaged over all sites of all files, and the clock rate has the same unit. The clock rate of file p is the clock rate of the run times μ_p; the model output reports it per file.

The normalization does not change the likelihood on branches at or above the floor. A branch below the floor, for example a zero-length branch, enters the propagation as μ_p * floor before the normalization and as (μ_p / m) * floor after it, so its term changes. The floor stays on the branch length, because single-file runs and v0 floor the branch length. The optimize loop therefore computes the likelihood again after each normalization. When the rates converge, m approaches 1 and the change vanishes.

**Rate step**: `fn optimize_gtr_rate()` ([packages/treetime/src/gtr/refinement.rs#L145-L253](../../packages/treetime/src/gtr/refinement.rs#L145-L253)) evaluates the likelihood at 0.01 * sqrt(μ), sqrt(μ) and 100 * sqrt(μ), and minimizes with a bracketed Brent search in sqrt(μ) when the middle point has the lowest cost. When it does not, the partition-rate step takes the end of the range with the lower cost (`RateRangeEdge::LowerCostEnd`), so a file whose likelihood keeps rising toward one end moves to that end. After the step, the caller runs the marginal update with the new rates, so the messages always belong to the current rates.

**Schedule**:

- **Optimize**: each iteration runs the branch-length step, the damping, the rate step and the normalization (`fn update_partition_rates()`, [packages/treetime/src/optimize/partition_rates.rs#L30-L48](../../packages/treetime/src/optimize/partition_rates.rs#L30-L48), called from [packages/treetime/src/optimize/run_loop.rs#L76-L83](../../packages/treetime/src/optimize/run_loop.rs#L76-L83)). The stop rule compares likelihoods after the normalization and needs no separate rate tolerance. The best iteration that the loop returns holds its branch lengths, rates and indel rates together
- **Timetree**: each of the two ML branch-length passes before the dating rounds ends with a rate step and a normalization ([packages/treetime/src/timetree/pre_loop.rs#L109-L116](../../packages/treetime/src/timetree/pre_loop.rs#L109-L116)). The dating rounds keep the rates fixed
- **Ancestral** keeps the input branch lengths, and **prune** counts parsimony substitutions only, so neither estimates rates

**Single file**: a run with one file has no rate step and no normalization, and `--partition-rates` has no effect. With `--model infer`, the model of a single-file optimize or timetree run has μ = 1 after the fit, as in v0 ([packages/legacy/treetime/treetime/treeanc.py#L1626-L1631](../../packages/legacy/treetime/treetime/treeanc.py#L1626-L1631)), in dense and sparse mode, so the branch lengths and the clock rate are expected substitutions per site. The rule keeps the branch lengths, as v0 does: optimize and timetree start the branch-length optimization from the input branch lengths. Ancestral reconstructs each file as its own partition on the input branch lengths, so it applies the single-file rule to every file, and a run with several files gives the models of separate runs. In optimize and timetree with several files, each file keeps its fitted μ as the start of the rate step; optimize also normalizes these fitted rates once before its first iteration ([packages/treetime/src/optimize/pipeline.rs#L94-L96](../../packages/treetime/src/optimize/pipeline.rs#L94-L96)). `enum InferredRate` ([packages/treetime/src/partition/inferred_rate.rs](../../packages/treetime/src/partition/inferred_rate.rs)) owns the rule for the first fit. The ancestral refinement of `--gtr-iterations` starts from that model and re-estimates μ with W and π ([kb/decisions/ancestral-iterative-gtr-refinement.md](ancestral-iterative-gtr-refinement.md)). Named models keep their own μ, and mugration keeps its fitted and optimized rate, as v0 does ([packages/legacy/treetime/treetime/wrappers.py#L787-L808](../../packages/legacy/treetime/treetime/wrappers.py#L787-L808)).

### Indels

- **One indel rate per file**: the indel rate λ_p of file p is the number of indel events of the file divided by the total branch length (`fn IndelPolicy::edge_indels()`, [packages/treetime/src/optimize/indel.rs#L24-L42](../../packages/treetime/src/optimize/indel.rs#L24-L42)). The model output reports it per file
- **One combined term**: all files share the branch lengths, so the sum of the per-file Poisson terms of a branch depends on the branch length as one Poisson term with the summed count and the summed rate Σ_p λ_p. TreeTime evaluates that one term, so every consumer of the indel term of [kb/decisions/optimize-indel-contribution-to-likelihood.md](optimize-indel-contribution-to-likelihood.md) works unchanged
- **File boundaries**: a deletion that crosses the boundary between two files counts as one event in each file. Terminal gap filling at the edge of a file can turn such gaps into unknown characters, so a split alignment gives the results of the unsplit alignment only when no gap touches the cut

## Why

- **Segmented genomes and genes**: influenza segments, or genes with different composition, evolve under different substitution models but share one tree. One file per partition is the input that segment builds already produce
- **Relative rates**: files evolve at different speeds. With shared branch lengths and one fixed μ per file, a slow file would pull the branch lengths down and a fast file would push them up. One relative rate per file with shared branch lengths is the default of RAxML-NG (`--brlen scaled`) and the `-p` option of IQ-TREE
- **Normalization**: the likelihood depends on μ_p and b only through the products μ_p * b, so the rates and the branch lengths have one free common scale. The site-weighted mean of 1 fixes that scale and gives branch lengths in substitutions per site averaged over all sites, the unit of a single-file run with μ = 1
- **Range end**: a file without any change prefers the lowest rate, so its optimum lies outside every range. Keeping the start rate in that case leaves a rate that the data contradict. When the middle point is not the lowest, the end with the lower cost costs at most as much as the start, so the step never lowers the likelihood
- **Indels**: one rate per file keeps the indel process of each file separate, as the substitution model is, and the combined term keeps one likelihood formula for all branch-length consumers

## Verification

v0 accepts one alignment and has no counterpart, so the evidence comes from tests:

- **Split equality**: one alignment and the same alignment cut into two column blocks give the same branch lengths, dates and pruned trees with `--model jc69` and `--partition-rates equal` when no gap touches the cut. Optimize and timetree are checked in dense and in sparse mode; prune always uses the sparse representation
- **Order**: the most likely results do not depend on the order of the files
- **Ancestral factorization**: for the most likely reconstruction, an ancestral run with two files equals two single-file runs, and its log-likelihood equals their sum. Posterior sampling (`--sample-from-profile` with `root` or `all`) draws from one random generator per run, file after file in the order of the `--alignment` list. A sampled result therefore depends on the order of the files, and every file after the first draws other values than in a separate run with the same seed
- **Scale**: the normalized rates and branch lengths do not depend on the scale of the starting rates, and one normalization keeps the likelihood on branches above the floor
- **Rate recovery**: with 20 replicates of three files of 3000 sites simulated on `data/zika/86/tree.nwk` under JC69 with the rates 0.3, 1.0 and 2.0, the mean estimate of each file lies within three standard errors of its normalized rate (3/11, 10/11, 20/11) (`packages/treetime/src/optimize/__tests__/test_pipeline_partition_rate_recovery.rs`)

## v0 behavior

v0 reads one alignment into `TreeAnc` with one GTR model ([packages/legacy/treetime/treetime/treeanc.py#L171-L184](../../packages/legacy/treetime/treetime/treeanc.py#L171-L184)). Users concatenate segments into one alignment, which gives every segment the same model and the same rate.

## Related

- [kb/decisions/partition-system-architecture.md](partition-system-architecture.md)
- [kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md)
- [kb/issues/M-partition-site-ranges-and-per-partition-models-unsupported.md](../issues/M-partition-site-ranges-and-per-partition-models-unsupported.md)
- [kb/issues/M-mugration-gtr-rate-restore-inconsistency.md](../issues/M-mugration-gtr-rate-restore-inconsistency.md)
- [kb/proposals/config-file-multi-partition.md](../proposals/config-file-multi-partition.md)
