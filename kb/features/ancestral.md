# Ancestral Reconstruction

## Fitch Parsimony (Sparse)

- [x] Backward pass (intersection/union of child state sets)
- [x] Forward pass (resolve ambiguities top-down)
- [x] Indel tracking (insertions, deletions)
- [x] Composition tracking (character counts)
- [x] Deterministic root state resolution ([intentional change](../decisions/ancestral-fitch-deterministic-root-state.md))
- [x] Forward cleanup drops stored full sequences from non-root internal nodes

## Marginal Reconstruction (Dense)

- [x] Backward pass (log-space message multiplication)
- [x] Forward pass (parent-child message propagation)
- [x] Root equilibrium frequency weighting
- [x] Log-space normalization (logsumexp)
- [x] expQt matrix propagation
- [x] `initialize_marginal()` bootstrap with dummy JC69
- [x] Pre-inference marginal pass to populate profiles
- [x] GTR selection after profiles exist
- [x] Second marginal pass after final GTR assignment
- [ ] Profile sampling (`sample_from_prof` - v0 supports stochastic sampling, v1 argmax only)

## Marginal Reconstruction (Sparse)

- [x] Backward pass (variable positions only)
- [x] Forward pass (variable positions only)
- [x] Fixed-position per-character profiles
- [x] Mutation composition from edge substitutions
- [x] Message combine with variability threshold (EPS=1e-4)
- [x] Sparse compression before GTR selection
- [x] Per-partition GTR assignment after dummy JC69 bootstrap
- [x] Single `update_marginal()` pass after final GTR assignment

## Joint Reconstruction

- [ ] Removed in v1 ([intentional change](../decisions/ancestral-joint-reconstruction-removed.md))
- [ ] Parsed and exposed in CLI but runtime path is `unimplemented!`

## GTR Model Inference

- [x] GTR inference from data (sparse and dense paths)
- [x] Uninformative root state filtering in dense path ([intentional change](../decisions/gtr-uninformative-root-state-filtering.md))
- [x] GTR bootstrapped through temporary JC69 model before replacement
- [ ] Iterative GTR inference (`infer_gtr_iterative()` in v0)

## CLI Options

- [x] `--reconstruct-tip-states` (overwrite ambiguous tips, controls leaf emission)
- [x] `--model` relevant only on marginal paths (parsimony bypasses GTR)
- [ ] `--keep-overhangs` (parsed but not wired; gap handling not implemented)
- [ ] `--zero-based` indexing (parsed but not wired)
- [ ] `--report-ambiguous` (parsed but not wired)
- [ ] `--seed` for reproducibility (parsed but not wired)
- [ ] `--gtr-params` custom GTR parameters (parsed but not wired)
- [ ] `--aa` (parsed but not wired)
- [ ] `--vcf-reference` (parsed, VCF reader not implemented)
- [ ] `--aln` legacy option (parsed but not wired)
- [ ] Tree inference from alignment (help text mentions it, not implemented)
