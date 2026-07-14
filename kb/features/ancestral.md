# Ancestral Reconstruction

## Fitch Parsimony (Sparse)

- [/] Backward pass (binary trees are covered; multifurcation recurrence can exceed the minimum score: [kb/issues/M-ancestral-fitch-polytomy-recurrence-not-minimum.md](../issues/M-ancestral-fitch-polytomy-recurrence-not-minimum.md))
- [x] Forward pass (resolve ambiguities top-down)
- [x] Indel tracking (insertions, deletions)
- [x] Composition tracking (character counts)
- [x] Deterministic root state resolution ([kb/decisions/ancestral-fitch-deterministic-root-state.md](../decisions/ancestral-fitch-deterministic-root-state.md))
- [x] Forward cleanup drops stored full sequences from non-root internal nodes
- [/] Parallel site classification ([kb/issues/N-ancestral-fitch-site-classification-parallel-regression.md](../issues/N-ancestral-fitch-site-classification-parallel-regression.md))

## Marginal Reconstruction (Dense)

- [/] Backward pass (impossible-factor multiplicity is not preserved by the cavity sentinel: [kb/issues/H-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md](../issues/H-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md))
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
- [/] Forward pass (variable positions only; reconstruction can remove state before a fallible lookup: [kb/issues/N-ancestral-sparse-remove-insert-pattern.md](../issues/N-ancestral-sparse-remove-insert-pattern.md))
- [x] Fixed-position per-character profiles
- [x] Mutation composition from edge substitutions
- [x] Message combine with variability threshold (EPS=1e-4)
- [x] Sparse compression before GTR selection
- [x] Per-partition GTR assignment after dummy JC69 bootstrap
- [x] Single `update_marginal()` pass after final GTR assignment
- [/] Fallible parallel passes can expose partial state before returning an error ([kb/issues/M-inference-fallible-parallel-passes-partially-commit.md](../issues/M-inference-fallible-parallel-passes-partially-commit.md))

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

- [x] `--reconstruct-tip-states` (overwrite ambiguous tips, controls leaf emission)
- [x] `--model` relevant only on marginal paths (parsimony bypasses GTR)
- [ ] `--keep-overhangs` (parsed but not wired; gap handling not implemented)
- [ ] `--zero-based` indexing (parsed but not wired; [kb/issues/M-core-mutation-representation-and-format-projection-inconsistent.md](../issues/M-core-mutation-representation-and-format-projection-inconsistent.md))
- [ ] `--report-ambiguous` (parsed but not wired)
- [x] `--seed` for reproducibility
- [ ] `--gtr-params` custom GTR parameters (parsed but not wired)
- [x] `--translations` per-CDS AA FASTA input (both `{cds}` and `%GENE` placeholders)
- [x] `--cdses` / `--genes` CDS names (derived from `--annotation` when omitted)
- [x] `--annotation` GFF3 annotation (CDS coordinates for node data; GenBank: [kb/issues/M-ancestral-genbank-annotation-unsupported.md](../issues/M-ancestral-genbank-annotation-unsupported.md))
- [x] `--aa-model` AA substitution model (infer default, jtt92 opt-in)
- [x] `--aa-root-sequence` per-CDS root/reference override
- [x] `--output-aa-sequences` per-CDS reconstructed AA FASTA output
- [ ] `--vcf-reference` (parsed, VCF reader not implemented)
- [ ] `--aln` legacy option (parsed but not wired)
- [ ] Tree inference from alignment (help text mentions it, not implemented)
