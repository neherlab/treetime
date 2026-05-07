# Mugration (Discrete Trait Reconstruction)

[Back to index](README.md)

Discrete trait ancestral reconstruction treats categorical metadata (locations, hosts, lineages) as characters evolving on the phylogeny. The term "mugration" (mutation + migration) reflects the model: transitions between discrete states are treated like substitutions using GTR-like machinery.

---

## Discrete Marginal Reconstruction

Felsenstein pruning (Felsenstein 1981) applied to discrete categorical traits rather than nucleotide sequences. The algorithm is identical to marginal ML ancestral reconstruction but operates on a discrete state alphabet (e.g., country names) with a GTR-like transition matrix.

v1: [`packages/treetime/src/commands/mugration/discrete_marginal.rs`](../../packages/treetime/src/commands/mugration/discrete_marginal.rs).
v0: [`packages/legacy/treetime/treetime/wrappers.py#L653-L811`](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811).

Key functions: `run_discrete_marginal()`, `attach_traits()`.

### Algorithm

**Initialization** (`attach_traits()` in [`packages/treetime/src/commands/mugration/discrete_marginal.rs#L96`](../../packages/treetime/src/commands/mugration/discrete_marginal.rs#L96)):

For each leaf node:

- **Known trait**: Create delta profile (one-hot encoding). `profile[i] = 1.0` where `i` is the observed state index, `0.0` elsewhere.
- **Missing trait ("?")**: Create uniform prior. `profile[i] = 1/n_states` for all states.

This follows Felsenstein's treatment of ambiguous data: unknown states receive equal probability across all possibilities, enabling marginalization during message passing.

**Backward pass** (postorder, leaves to root) (`process_node_backward()` in [`packages/treetime/src/representation/partition/discrete.rs#L60`](../../packages/treetime/src/representation/partition/discrete.rs#L60)):

For each node in postorder:

1. If leaf: use initialized profile directly
2. If internal: combine child messages via element-wise product in log space
   ```
   log_profile = sum_c log(msg_from_child_c)
   profile = normalize(exp(log_profile))
   ```
3. If root: weight by equilibrium frequencies `pi`
   ```
   profile = normalize(profile * pi)
   ```
4. If not root: propagate message to parent via transition matrix
   ```
   msg_from_child[j] = sum_i profile[i] * P(i->j|t)
   ```
   where `P = exp(Q*t)` is the transition probability matrix for branch length `t`.

**Forward pass** (preorder, root to leaves) (`process_node_forward()` in [`packages/treetime/src/representation/partition/discrete.rs#L150`](../../packages/treetime/src/representation/partition/discrete.rs#L150)):

For each node in preorder:

1. If root: profile already computed
2. If non-root: combine incoming message from parent with message to parent
   ```
   profile = normalize(msg_to_parent * propagated_msg_from_parent)
   ```
   where `propagated_msg_from_parent = exp(Q*t) . msg_to_child`
3. Compute outgoing `msg_to_child` for each child edge:
   ```
   msg_to_child = normalize(profile / msg_from_child)
   ```

**Trait assignment** (`get_reconstructed_trait()` in [`packages/treetime/src/representation/partition/discrete.rs#L50`](../../packages/treetime/src/representation/partition/discrete.rs#L50)):

After forward-backward passes, each node has a posterior probability distribution over states. The assigned trait is `argmax(profile)`.

### Missing Data Handling

Leaves with missing traits (`"?"` in metadata) are valid reconstruction targets:

- Uniform prior allows tree context to inform inference
- After forward-backward, the posterior reflects phylogenetic signal from neighboring nodes
- Assigned trait is the most likely state given tree structure

This is analogous to `--reconstruct-tip-states` in the ancestral command but applied to discrete traits rather than nucleotide ambiguities.

Supporting references for missing data treatment:

- Felsenstein (2003), p. 255: ambiguous states receive likelihood 1.0 for all possibilities
- PastML (Ishikawa et al. 2019, Mol Biol Evol): "unknown states can be omitted and will be estimated during analysis"
- BEAST UTM model (Vaiente & Scotch 2020): extends to prior probability distributions over uncertain tip states

---

## GTR Model Construction

Constructs a GTR-like transition model for discrete traits.

v1: [`packages/treetime/src/commands/mugration/run.rs#L223`](../../packages/treetime/src/commands/mugration/run.rs#L223).

**v1 implementation**:

- Equilibrium frequencies `pi`: uniform (1/n_states) or from weights file, with pseudo-count smoothing
- Exchangeability matrix `W`: uniform (all transitions equally likely)
- Initial forward-backward pass with uniform model, then iterative GTR refinement via `refine_gtr_iterative()` ([`packages/treetime/src/commands/mugration/gtr_refinement.rs#L28`](../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L28))

**v0 implementation** (see [Iterative GTR for Discrete Traits](unimplemented.md#iterative-gtr-for-discrete-traits-ported)):

- Initial reconstruction with uniform model
- 5 iterations of `infer_gtr()` + `optimize_gtr_rate()` re-estimation
- Final reconstruction with refined model

Both v0 and v1 perform iterative GTR refinement. v1's implementation includes two intentional improvements (pseudo-count smoothing on initial pi, root state uniform-threshold filtering) that shift posterior probabilities at ambiguous internal nodes. Remaining parity differences tracked in [Mugration golden master parity with v0](../issues/M-mugration-iterative-gtr.md).

---

## Confidence Profiles

After forward-backward, `get_confidence(node_key)` returns the full posterior distribution `profile[n_states]` at [`packages/treetime/src/representation/partition/discrete.rs#L56`](../../packages/treetime/src/representation/partition/discrete.rs#L56). This enables uncertainty quantification for trait assignments and identification of ambiguous nodes (flat profiles).

---

## Output Files

| File                   | Content                                        | v1 Status |
| ---------------------- | ---------------------------------------------- | --------- |
| `traits.csv`           | Per-node trait assignments (all nodes)         | Complete  |
| `annotated_tree.nexus` | Newick tree with trait annotations in comments | Complete  |
| `annotated_tree.nwk`   | Newick tree with NHX-style trait annotations   | Complete  |
| `gtr.json`             | GTR model parameters                           | Complete  |
| `confidence.csv`       | Per-node posterior distributions (optional)    | Complete  |

---

## Test Coverage

| Test Type     | Location                                                                                                          | Coverage                               |
| ------------- | ----------------------------------------------------------------------------------------------------------------- | -------------------------------------- |
| Golden master | [`__tests__/test_gm_mugration.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_gm_mugration.rs) | Internal nodes only (filtered)         |
| Command       | [`__tests__/test_run.rs`](../../packages/treetime/src/commands/mugration/__tests__/test_run.rs)                   | Output file existence, basic structure |
| Partition     | [`representation/partition/discrete.rs`](../../packages/treetime/src/representation/partition/discrete.rs)        | Unit tests for node data, profiles     |

### Test Gaps

1. **Leaf nodes filtered from golden master comparison**: Test compares only `NODE_*` prefixed names
2. **Missing data reconstruction not tested**: No cases with `"?"` trait values
3. **Confidence profiles not validated**: `MugrationOutput.confidence` captured but not asserted
4. **Degenerate cases not covered**: Single state, all missing, etc.

---

## Known Issues

- [Mugration golden master parity with v0](../issues/M-mugration-iterative-gtr.md) - iterative GTR implemented but two intentional v1 improvements cause argmax divergence at ambiguous internal nodes for 5/7 test datasets

---

## References

### Primary

Sagulenko, Puller & Neher (2018). "TreeTime: Maximum-likelihood phylodynamic analysis." Virus Evolution, 4(1):vex042. doi:10.1093/ve/vex042

### Ancestral Reconstruction

- Felsenstein (1981). "Evolutionary trees from DNA sequences: a maximum likelihood approach." J Mol Evol, 17(6):368-376. doi:10.1007/BF01734359
- Felsenstein (2003). "Inferring Phylogenies." Sinauer Associates, p. 255.

### Phylogeography

- Lemey, Rambaut, Drummond & Suchard (2009). "Bayesian phylogeography finds its roots." PLoS Computational Biology, 5(9):e1000520.
- Edwards et al. (2011). "Ancient hybridization and an Irish origin for the modern polar bear matriline." Current Biology, 21(15):1251-1258.

### Missing Data

- Ishikawa et al. (2019). "A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios." Mol Biol Evol, 36(9):2069-2085.
- Vaiente, M.A. & Scotch, M. (2020). "Going back to the roots: Evaluating Bayesian phylogeographic models with discrete trait uncertainty." Infection, Genetics and Evolution, 85:104501. doi:10.1016/j.meegid.2020.104501
