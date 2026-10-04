# Mugration (Discrete Trait Reconstruction)

Discrete trait ancestral reconstruction treats categorical metadata (locations, hosts, lineages) as characters evolving on the phylogeny. The term "mugration" (mutation + migration) reflects the model: transitions between discrete states are treated like substitutions using GTR-like machinery.

---

## Discrete Marginal Reconstruction

Felsenstein pruning (<a id="cite-1"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[1](#ref-1)]) applied to discrete categorical traits rather than nucleotide sequences. The algorithm is identical to marginal ML ancestral reconstruction but operates on a discrete state alphabet (e.g., country names) with a GTR-like transition matrix.

v1: shared marginal passes in [`packages/treetime/src/partition/marginal/shared/`](../../packages/treetime/src/partition/marginal/shared/) with the discrete partition in [`packages/treetime/src/partition/marginal/discrete/partition.rs`](../../packages/treetime/src/partition/marginal/discrete/partition.rs).
v0: [`packages/legacy/treetime/treetime/wrappers.py#L653-L811`](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811).

Key functions: `MarginalPasses::marginal_update()`, `PartitionMarginalDiscrete::new()`.

### Algorithm

Initialization (`PartitionMarginalDiscrete::new()` in [`packages/treetime/src/partition/marginal/discrete/partition.rs#L27`](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L27)) stores one observation per leaf, and every backward pass starts from it:

For each leaf node:

- Known trait: Create delta profile (one-hot encoding). `profile[i] = 1.0` where `i` is the observed state index, `0.0` elsewhere.
- Missing trait ("?"): `profile[i] = 1.0` for all states, as v0 (`mugration_GTR.profile_map[missing_char] = np.ones(n_states)`, [`packages/legacy/treetime/treetime/wrappers.py#L767`](../../packages/legacy/treetime/treetime/wrappers.py#L767)).

This follows Felsenstein's treatment of ambiguous data: the likelihood of an unknown state is 1 for every possible state, so the leaf contributes a factor of 1 and the tree likelihood marginalizes over it. A profile of `1/n_states` gives the same posteriors but lowers the log-likelihood by `ln(n_states)` per missing leaf.

Backward pass (postorder, leaves to root) (`fn indexed_node_backward` [packages/treetime/src/partition/marginal/shared/pass.rs#L83](../../packages/treetime/src/partition/marginal/shared/pass.rs#L83)):

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

Forward pass (preorder, root to leaves) (`fn indexed_node_forward` [packages/treetime/src/partition/marginal/shared/pass.rs#L225](../../packages/treetime/src/partition/marginal/shared/pass.rs#L225)):

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

Trait assignment (`fn get_reconstructed_trait` [packages/treetime/src/partition/marginal/discrete/partition.rs#L67](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L67)):

After forward-backward passes, each node has a posterior probability distribution over states. The assigned trait is `argmax(profile)`.

### Missing Data Handling

Leaves with missing traits (`"?"` in metadata) are valid reconstruction targets:

- Uniform prior allows tree context to inform inference
- After forward-backward, the posterior reflects phylogenetic signal from neighboring nodes
- Assigned trait is the most likely state given tree structure

This is analogous to `--reconstruct-tip-states` in the ancestral command but applied to discrete traits rather than nucleotide ambiguities.

Supporting references for missing data treatment:

- <a id="cite-2"></a>[Felsenstein 2003](https://doi.org/10.1007/978-0-387-21337-7) [[2](#ref-2)], p. 255: ambiguous states receive likelihood 1.0 for all possibilities
- PastML (<a id="cite-3"></a>[Ishikawa et al. 2019](https://doi.org/10.1093/molbev/msz131) [[3](#ref-3)]): "unknown states can be omitted and will be estimated during analysis"
- BEAST UTM model (<a id="cite-4"></a>[Vaiente and Scotch 2020](https://doi.org/10.1016/j.meegid.2020.104501) [[4](#ref-4)]): extends to prior probability distributions over uncertain tip states

---

## Architecture

The discrete partition (`struct PartitionMarginalDiscrete`) implements the shared `trait MarginalPasses` [packages/treetime/src/partition/marginal/shared/update.rs#L11](../../packages/treetime/src/partition/marginal/shared/update.rs#L11). Its passes call the shared `fn indexed_backward` and `fn indexed_forward` with `IndexedKind::Discrete` [packages/treetime/src/partition/marginal/discrete/partition.rs#L87-L91](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L87-L91). Discrete traits are `Array2(1, n_states)` profiles stored in the same `DenseNodeState` type as dense sequence inference. A discrete leaf reads its observed profile in `fn discrete_leaf_backward` [packages/treetime/src/partition/marginal/shared/pass.rs#L160](../../packages/treetime/src/partition/marginal/shared/pass.rs#L160). The forward pass skips the indel and sequence step for the discrete kind [packages/treetime/src/partition/marginal/shared/pass.rs#L263-L268](../../packages/treetime/src/partition/marginal/shared/pass.rs#L263-L268).

GTR refinement lives in [`packages/treetime/src/gtr/refinement.rs`](../../packages/treetime/src/gtr/refinement.rs) as a domain-level algorithm generic over `MarginalPasses`, called from the mugration pipeline.

---

## GTR Model Construction

Constructs a GTR-like transition model for discrete traits.

v1: `fn run` [packages/treetime/src/mugration/pipeline.rs#L84-L97](../../packages/treetime/src/mugration/pipeline.rs#L84-L97).

v1 implementation:

- Equilibrium frequencies `pi`: uniform (1/n_states) or from weights file, with pseudo-count smoothing
- Exchangeability matrix `W`: uniform (all transitions equally likely)
- Initial forward-backward pass with uniform model, then iterative GTR refinement via `fn refine_gtr_model_and_rate` [packages/treetime/src/gtr/refinement.rs#L50](../../packages/treetime/src/gtr/refinement.rs#L50)

v0 implementation (see [Iterative GTR for Discrete Traits](unimplemented.md#iterative-gtr-for-discrete-traits-ported)):

- Initial reconstruction with uniform model
- 5 iterations of `infer_gtr()` + `optimize_gtr_rate()` re-estimation
- Final reconstruction with refined model

Both v0 and v1 perform iterative GTR refinement. v1's implementation includes two intentional improvements (pseudo-count smoothing on initial pi, root state uniform-threshold filtering) that shift posterior probabilities at ambiguous internal nodes. Remaining parity differences tracked in [Mugration golden master parity with v0](../issues/M-mugration-iterative-gtr.md).

---

## Confidence Profiles

After forward-backward, `fn get_confidence` [packages/treetime/src/partition/marginal/discrete/partition.rs#L78](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L78) returns the full posterior distribution `profile[n_states]` of a node. This enables uncertainty quantification for trait assignments and identification of ambiguous nodes (flat profiles).

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

| Test Type     | Location                                                                                                           | Coverage                                           |
| ------------- | ------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------- |
| Golden master | [`__tests__/test_gm_mugration.rs`](../../packages/treetime/src/mugration/__tests__/test_gm_mugration.rs)           | Internal nodes only (filtered)                     |
| Command       | [`__tests__/test_run.rs`](../../packages/treetime/src/mugration/__tests__/test_run.rs)                             | Output file existence, basic structure             |
| Partition     | [`__tests__/test_discrete_marginal.rs`](../../packages/treetime/src/mugration/__tests__/test_discrete_marginal.rs) | Unit tests for attach, backward/forward, full pass |
| Brent         | inline tests in [`gtr/brent_bracketed.rs`](../../packages/treetime/src/gtr/brent_bracketed.rs)                     | Pure math optimizer tests                          |

### Test Gaps

1. Leaf nodes filtered from golden master comparison: Test compares only `NODE_*` prefixed names
2. Missing data reconstruction not tested: No cases with `"?"` trait values
3. Confidence profiles not validated: `MugrationOutput.confidence` captured but not asserted
4. Degenerate cases not covered: Single state, all missing, etc.
5. `count_transitions()` has no direct unit test ([N-mugration-count-transitions-untested](../issues/N-mugration-count-transitions-untested.md))

---

## Known Issues

- [Mugration golden master parity with v0](../issues/M-mugration-iterative-gtr.md) - iterative GTR implemented but two intentional v1 improvements cause argmax divergence at ambiguous internal nodes for 5/7 test datasets
- [Mugration count_transitions lacks direct unit test](../issues/N-mugration-count-transitions-untested.md)

---

## References

- <a id="ref-0"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042
- <a id="ref-1"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-1)
- <a id="ref-2"></a>Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. [↩](#cite-2)
- <a id="ref-3"></a>Ishikawa, Sohta A., Anna Zhukova, Wataru Iwasaki, and Olivier Gascuel. 2019. "A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios." _Molecular Biology and Evolution_ 36(9):2069-2085. https://doi.org/10.1093/molbev/msz131 [↩](#cite-3)
- <a id="ref-4"></a>Vaiente, Mara A., and Matthew Scotch. 2020. "Going Back to the Roots: Evaluating Bayesian Phylogeographic Models with Discrete Trait Uncertainty." _Infection, Genetics and Evolution_ 85:104501. https://doi.org/10.1016/j.meegid.2020.104501 [↩](#cite-4)
- <a id="ref-5"></a>Lemey, Philippe, Andrew Rambaut, Alexei J. Drummond, and Marc A. Suchard. 2009. "Bayesian Phylogeography Finds Its Roots." _PLoS Computational Biology_ 5(9):e1000520. https://doi.org/10.1371/journal.pcbi.1000520
- <a id="ref-6"></a>Edwards, Ceiridwen J., Marc A. Suchard, Philippe Lemey, et al. 2011. "Ancient Hybridization and an Irish Origin for the Modern Polar Bear Matriline." _Current Biology_ 21(15):1251-1258. https://doi.org/10.1016/j.cub.2011.05.058
