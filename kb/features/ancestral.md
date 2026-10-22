# Ancestral Reconstruction

## Fitch Parsimony (Sparse)

- [x] Backward pass (multifurcations use the plurality recurrence, a deliberate divergence from v0: [kb/decisions/ancestral-fitch-plurality-on-multifurcations.md](../decisions/ancestral-fitch-plurality-on-multifurcations.md))
- [x] Forward pass (resolve ambiguities top-down)
- [x] Indel tracking (insertions, deletions)
- [x] Composition tracking (character counts)
- [x] Deterministic root state resolution ([kb/decisions/ancestral-fitch-deterministic-root-state.md](../decisions/ancestral-fitch-deterministic-root-state.md))
- [x] Forward cleanup drops stored full sequences from non-root internal nodes
- [/] Fitch preprocessing parallelism (the former serial site index is gone and leaf construction is parallel; site-level work and complete-command scaling remain unverified: [kb/issues/N-ancestral-fitch-site-classification-parallel-scaling-unverified.md](../issues/N-ancestral-fitch-site-classification-parallel-scaling-unverified.md))

## Marginal Reconstruction (Dense)

- [/] Backward pass (impossible-factor multiplicity is not preserved by the cavity sentinel: [kb/issues/M-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md](../issues/M-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md))
- [x] Forward pass (parent-child message propagation)
- [x] Root equilibrium frequency weighting
- [x] Log-space normalization (logsumexp)
- [x] expQt matrix propagation
- [x] `initialize_marginal()` bootstrap with dummy JC69
- [x] Pre-inference marginal pass to populate profiles
- [x] GTR selection after profiles exist
- [x] Second marginal pass after final GTR assignment
- [x] Profile sampling (`--sample-from-profile` with argmax/root/all modes, seeded RNG)

## Marginal Reconstruction (Sparse)

- [x] Backward pass (variable positions only)
- [x] Forward pass (variable positions only)
- [x] Fixed-position per-character profiles
- [x] Mutation composition from edge substitutions
- [x] Message combine with variability threshold (EPS=1e-4)
- [x] Sparse compression before GTR selection
- [x] Per-partition GTR assignment after dummy JC69 bootstrap
- [x] Single `update_marginal()` pass after final GTR assignment

## Several Alignment Files

- [x] One partition per `--alignment` file, each with its own model and reconstruction ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))
- [x] For the most likely reconstruction, a run with several files equals separate single-file runs, and the log-likelihood is their sum. Posterior sampling (`--sample-from-profile` with `root` or `all`) draws from one random generator per run, file after file in the order of the `--alignment` list. A sampled result therefore depends on the order of the files, and every file after the first draws other values than in a separate run with the same seed
- [x] Each file keeps its own short-branch floor
- [x] Missing-tip check per file (stop above one third unless `--allow-missing-sequences`)
- [x] Per-file sequence outputs of several alignment files (FASTA, Auspice JSON, UShER MAT, Nexus, Newick with `--newick-annotations beast` or `nhx`, ancestral augur node data): `<command>.<n>.<ext>` under `--output-all`, or a `{alignment}` placeholder in the per-file flag; a per-file flag without the placeholder stops with an error
- [x] `--translations` and `--annotation` stop with an error on several files

## Per-CDS Amino-Acid Reconstruction

- [x] Multi-partition marginal traversal on shared graph (one parse, key-join)
- [x] Stop-inclusive AA alphabet (22-state `aa`) for default infer model
- [x] Empirical model out-of-alphabet handling (map to unknown, warn) ([kb/decisions/ancestral-aa-empirical-model-out-of-alphabet-to-unknown.md](../decisions/ancestral-aa-empirical-model-out-of-alphabet-to-unknown.md))
- [x] GFF3 CDS parsing with nextclade name-resolution priority
- [x] Minus-strand compound CDS segments ordered 5'-to-3'
- [x] CDS nucleotide length divisible-by-3 validation
- [x] AA-to-annotation length invariant (`3 * aa_len == cds_nt_len`)
- [x] Augur node-data JSON: `aa_muts`, `aa_sequences` (root-only), per-CDS annotations, per-CDS reference
- [x] Missing-tip handling shared with nuc (all-unknown synthesis, >1/3 threshold)
- [ ] Golden master validation vs augur (blocked: augur not installed in dev container)

## Joint Reconstruction

- [ ] Removed in v1 ([kb/decisions/ancestral-joint-reconstruction-removed.md](../decisions/ancestral-joint-reconstruction-removed.md))
- [ ] Parsed and exposed in CLI but runtime path is `unimplemented!`

## GTR Model Inference

- [x] GTR inference from data (sparse and dense paths)
- [x] Uninformative root state filtering in dense path ([kb/decisions/gtr-uninformative-root-state-filtering.md](../decisions/gtr-uninformative-root-state-filtering.md))
- [x] GTR bootstrapped through temporary JC69 model before replacement
- [ ] Iterative GTR inference (`infer_gtr_iterative()` in v0)

## CLI Options

- [x] `--include-leaves` (emit reconstructed tip sequences)
- [x] `--impute-missing-data` (resolve ambiguous/unknown tips to the marginal argmax, leaving gaps; no-op under parsimony). See [kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md)
- [x] `--reconstruct-tip-states` (v0-compatible alias enabling both `--include-leaves` and `--impute-missing-data`)
- [x] `--model` and `--model-params` apply only with `--method-anc marginal`; with parsimony they stop the run with an argument error (parsimony uses no model)
- [x] `--zero-based` (mutation positions in the comments of Newick and Nexus trees)
- [x] `--report-ambiguous` (keep mutations from or to `N`/`X` in written mutation lists; omitted by default, as in v0; see [kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md))
- [x] `--seed` for reproducibility; a generated seed is logged when sampling runs without one (see [kb/decisions/cli-seed-on-commands-with-random-steps.md](../decisions/cli-seed-on-commands-with-random-steps.md))
- [x] `--model-params` (`key=value` parameters of a named model, as v0 `--gtr-params`; an unknown key or a malformed item is an error)
- [x] `--translations` per-CDS AA FASTA input (both `{cds}` and `%GENE` placeholders)
- [x] `--cds-selection` CDS names (derived from `--annotation` when omitted)
- [x] `--annotation` GFF3 annotation (CDS coordinates for node data; GenBank: [kb/issues/M-ancestral-genbank-annotation-unsupported.md](../issues/M-ancestral-genbank-annotation-unsupported.md))
- [x] `--aa-model` AA substitution model (infer default, jtt92 opt-in)
- [x] `--aa-root-sequence` per-CDS root/reference override
- [x] `--output-aa-sequences` per-CDS reconstructed AA FASTA output
- [x] Optional ambiguity report as IUPAC codes in the reconstructed FASTA outputs and the node sequences of augur node data (`--ambiguity-threshold`, `--ambiguity-outputs`; new in v1, [kb/decisions/ambiguity-report-of-plausible-states.md](../decisions/ambiguity-report-of-plausible-states.md))
- [ ] VCF input with `--vcf-reference` ([kb/issues/M-io-vcf-input-output-unimplemented.md](../issues/M-io-vcf-input-output-unimplemented.md))
- [ ] `--aln` legacy option (parsed but not wired)
- [ ] Tree inference from alignment (help text mentions it, not implemented)
