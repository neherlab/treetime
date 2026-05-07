# GTR Substitution Models

## Implemented Models

- [x] Named nucleotide models: JC69, K80, F81, HKY85, T92, TN93
- [x] Amino acid model: JTT92

## Model Inference

- [x] Model inference from data (sparse and dense paths)
- [x] Eigendecomposition (symmetric, caching)
- [x] expQt (matrix exponential for transition probabilities)
- [x] Equilibrium frequencies (pi vector, normalized)
- [x] Exchangeability matrix (W, symmetric, normalized by avg_transition)
- [x] Rate scaling (mu parameter)
- [x] JSON output (GtrOutput struct with model type, parameters)

## Not Implemented

- [ ] Custom GTR from file (`--custom-gtr` / `GTR.from_file()` in v0)
- [ ] Random GTR generation (`GTR.random()` in v0, used for testing)
- [ ] GTR save to file (`save_to_npz()` in v0)
- [ ] Site-specific models (multi-dimensional pi/W in v0)
- [ ] Branch length optimization at GTR level (v0 `optimal_t()` - moved to partition layer in v1)
- [ ] Sequence probability at GTR level (v0 `prob_t()` - moved to partition layer in v1)
- [ ] State pair compression (v0 `state_pair()` - moved to partition layer in v1)
