# Mugration (Discrete Trait Reconstruction)

Discrete trait ancestral reconstruction treats categorical metadata (locations, hosts, lineages) as characters evolving on the phylogeny. The term "mugration" (mutation + migration) reflects the model: transitions between discrete states are treated like substitutions using GTR-like machinery.

---

## Discrete Marginal Reconstruction

Felsenstein pruning (<a id="cite-1"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[1](#ref-1)]) applied to discrete categorical traits rather than nucleotide sequences. The algorithm is identical to marginal ML ancestral reconstruction but operates on a discrete state alphabet (e.g., country names) with a GTR-like transition matrix.

v1: shared marginal passes in [`packages/treetime/src/partition/marginal/shared/`](../../packages/treetime/src/partition/marginal/shared/) with the discrete partition in [`packages/treetime/src/partition/marginal/discrete/partition.rs`](../../packages/treetime/src/partition/marginal/discrete/partition.rs), reached through the marginal fit of ancestral reconstruction (see [Architecture](#architecture)).
v0: [`packages/legacy/treetime/treetime/wrappers.py#L653-L811`](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811).

Key functions: `fit_marginal()`, `refine_gtr()`, `PartitionMarginalDiscrete::new()`, `PartitionMarginalDiscrete::state_posteriors()`.

### Algorithm

The trait states are a `struct StateAlphabet` of names ([`packages/treetime/src/alphabet/state_alphabet.rs#L10`](../../packages/treetime/src/alphabet/state_alphabet.rs#L10)): sorted, unique, with at least 2 states. The states are names, not ASCII characters, so the number of states has no limit other than memory, and non-ASCII names (`Côte d'Ivoire`, `北京`) are valid. The GTR `pi` and `W`, the profiles, and the outputs use the order of the state alphabet. The missing-data value (`--missing-data`, default `?`) is a metadata convention of the mugration front and is not a state.

`fn leaf_observations()` ([`packages/treetime/src/mugration/leaf_observations.rs#L13`](../../packages/treetime/src/mugration/leaf_observations.rs#L13)) gives each leaf an `enum StateObservation { Known(usize), Unknown }` ([`packages/treetime/src/partition/marginal/discrete/partition.rs#L84`](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L84)). Initialization (`PartitionMarginalDiscrete::new()` in [`packages/treetime/src/partition/marginal/discrete/partition.rs#L21`](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L21), profiles from `fn leaf_profile()` at [#L89-L98](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L89-L98)) stores one profile per leaf, and every backward pass starts from it:

For each leaf node:

- Known trait: Create delta profile (one-hot encoding). `profile[i] = 1.0` where `i` is the observed state index, `0.0` elsewhere.
- Unknown trait (the missing-data value, or no metadata row): `profile[i] = 1.0` for all states, as v0 (`mugration_GTR.profile_map[missing_char] = np.ones(n_states)`, [`packages/legacy/treetime/treetime/wrappers.py#L767`](../../packages/legacy/treetime/treetime/wrappers.py#L767)).

This follows Felsenstein's treatment of ambiguous data: the likelihood of an unknown state is 1 for every possible state, so the leaf contributes a factor of 1 and the tree likelihood marginalizes over it. A profile of `1/n_states` gives the same posteriors but lowers the log-likelihood by `ln(n_states)` per missing leaf.

Backward pass (postorder, leaves to root) (`fn IndexedBackwardPass::node` [packages/treetime/src/partition/marginal/shared/pass.rs#L150-L213](../../packages/treetime/src/partition/marginal/shared/pass.rs#L150-L213)):

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

Forward pass (preorder, root to leaves) (`fn IndexedForwardPass::node` [packages/treetime/src/partition/marginal/shared/pass.rs#L311-L353](../../packages/treetime/src/partition/marginal/shared/pass.rs#L311-L353)):

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

Trait assignment (`fn PartitionMarginalDiscrete::state_posteriors()` [packages/treetime/src/partition/marginal/discrete/partition.rs#L48-L68](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L48-L68)):

After the final marginal update, each node has a posterior probability distribution over states. `StatePosterior::new()` ([packages/treetime/src/partition/marginal/discrete/state_posterior.rs#L20-L32](../../packages/treetime/src/partition/marginal/discrete/state_posterior.rs#L20-L32)) computes one `struct StatePosterior { state, profile, entropy }` per node: the assigned trait is `argmax(profile)` (the first maximum, `argmax_first()`), and the entropy is `-Σ p ln(p + 1e-12)`. Every node has a posterior, leaves with an unknown trait included (v0 `reconstruct_tip_states=True`).

### Missing Data Handling

Leaves with unknown traits are valid reconstruction targets. A leaf is unknown in two cases:

- The metadata row holds the missing-data value: explicit missing data, no warning
- The leaf has no metadata row: the run continues, with one warning per leaf and one warning with their count from `fn warn_leaves_without_data()` ([packages/treetime/src/ancestral/attach.rs#L67-L89](../../packages/treetime/src/ancestral/attach.rs#L67-L89)), the warning owner that the sequence commands use for leaves without a sequence. As in v0 ([packages/legacy/treetime/treetime/wrappers.py#L777-L782](../../packages/legacy/treetime/treetime/wrappers.py#L777-L782)), no limit applies to the share of such leaves; the one-third limit of the sequence commands is a separate step of `fn complete_alignment_for_leaves()`

Metadata rows are paired with leaves by name in `fn pair_trait_rows()` ([packages/treetime/src/mugration/trait_pairing.rs#L6](../../packages/treetime/src/mugration/trait_pairing.rs#L6)), which returns a `struct TraitPairing` with the trait of each leaf, the observed values (unmatched rows included), the repeated names, and the unmatched rows.

For an unknown leaf:

- The all-ones profile lets the tree context inform the inference
- After forward-backward, the posterior reflects phylogenetic signal from neighboring nodes
- Assigned trait is the most likely state given tree structure

This is analogous to `--reconstruct-tip-states` in the ancestral command but applied to discrete traits rather than nucleotide ambiguities.

Supporting references for missing data treatment:

- <a id="cite-2"></a>[Felsenstein 2003](https://doi.org/10.1007/978-0-387-21337-7) [[2](#ref-2)], p. 255: ambiguous states receive likelihood 1.0 for all possibilities
- PastML (<a id="cite-3"></a>[Ishikawa et al. 2019](https://doi.org/10.1093/molbev/msz131) [[3](#ref-3)]): "unknown states can be omitted and will be estimated during analysis"
- BEAST UTM model (<a id="cite-4"></a>[Vaiente and Scotch 2020](https://doi.org/10.1016/j.meegid.2020.104501) [[4](#ref-4)]): extends to prior probability distributions over uncertain tip states

---

## Architecture

The discrete partition (`struct PartitionMarginalDiscrete`) implements the shared `trait MarginalPasses` [packages/treetime/src/partition/marginal/shared/update.rs#L10](../../packages/treetime/src/partition/marginal/shared/update.rs#L10). Its passes call the shared `fn indexed_backward` and `fn indexed_forward` with `IndexedKind::Discrete` [packages/treetime/src/partition/marginal/discrete/partition.rs#L76-L80](../../packages/treetime/src/partition/marginal/discrete/partition.rs#L76-L80). Discrete traits are `Array2(1, n_states)` profiles stored in the same `DenseNodeState` type as dense sequence inference. A discrete leaf reads its observed profile in `fn discrete_leaf_backward` [packages/treetime/src/partition/marginal/shared/pass.rs#L233-L240](../../packages/treetime/src/partition/marginal/shared/pass.rs#L233-L240). The forward pass skips the indel and sequence step for the discrete kind [packages/treetime/src/partition/marginal/shared/pass.rs#L344-L349](../../packages/treetime/src/partition/marginal/shared/pass.rs#L344-L349).

Mugration is an ancestral reconstruction of one discrete trait. `fn run()` ([packages/treetime/src/mugration/pipeline.rs#L26-L145](../../packages/treetime/src/mugration/pipeline.rs#L26-L145)) builds only what is specific to traits: the state alphabet from the observed values and the weight keys, the weight coverage check, `pi`, the start GTR, the leaf observations, and the partition. It then calls the marginal fit of ancestral reconstruction, `fn fit_marginal()` ([packages/treetime/src/ancestral/plan.rs#L292-L309](../../packages/treetime/src/ancestral/plan.rs#L292-L309)), which the dense and sparse sequence paths also call: a marginal update under the start model, then the optional GTR refinement. Last, it reads the state posteriors. `struct MugrationOutput { graph, gtr, states, posteriors }` carries the result to the writers.

The ancestral command types (`enum MarginalReconstruction`, `enum AncestralPartition`, `struct PartitionResult`) stay sequence-only: the discrete partition reaches the shared code through the generic `fit_marginal()`, as the amino-acid path reaches the ancestral plan.

GTR refinement lives in `fn refine_gtr()` ([packages/treetime/src/gtr/refinement.rs#L18-L95](../../packages/treetime/src/gtr/refinement.rs#L18-L95)), one routine generic over `MarginalPasses` for ancestral `--gtr-iterations` and mugration. A `struct GtrRefinement { iterations, pc, fixed_pi, rate }` selects the run; `enum RefinementRate` is `Keep` (ancestral) or `Optimize { scale }` (mugration).

The writers of every output (traits CSV, confidence CSV, augur node data, Auspice JSON, Newick and Nexus comments) read the state posteriors through `struct TreeTraits { attribute, states, posteriors }` ([packages/app-output/src/annotated_graph.rs#L86-L91](../../packages/app-output/src/annotated_graph.rs#L86-L91)) and compute no argmax or entropy. The branch transition label (`parent → child`) has one owner, `fn trait_change()` ([packages/app-output/src/trait_profile.rs#L16-L27](../../packages/app-output/src/trait_profile.rs#L16-L27)), and the confidence filter (`> 0.001`) has one owner, `fn build_confidence_map()` ([packages/app-output/src/trait_profile.rs#L9-L14](../../packages/app-output/src/trait_profile.rs#L9-L14)).

---

## GTR Model Construction

Constructs a GTR-like transition model for discrete traits.

v1: `fn run` [packages/treetime/src/mugration/pipeline.rs#L83-L130](../../packages/treetime/src/mugration/pipeline.rs#L83-L130).

Start model:

- Equilibrium frequencies `pi`: uniform (1/n_states) or from the weights file (fixed during refinement); optional pseudo-count smoothing with `--smooth-initial-pi`
- Exchangeability matrix `W`: uniform (all transitions equally likely), rate `Mu::FromMatrix(1.0)`

Refinement spec: `GtrRefinement { iterations: --iterations, pc: --pc (default 1.0), fixed_pi: pi of the weights, rate: Optimize { scale: --sampling-bias-correction } }`.

Refinement order of `fn refine_gtr()` with `RefinementRate::Optimize`, the order of v0's `reconstruct_discrete_traits()` ([packages/legacy/treetime/treetime/wrappers.py#L785-L811](../../packages/legacy/treetime/treetime/wrappers.py#L785-L811)):

1. A marginal update under the start model (in `fit_marginal()`), and the first fit from it
2. A full marginal update under the first fitted model ([packages/treetime/src/gtr/refinement.rs#L55-L62](../../packages/treetime/src/gtr/refinement.rs#L55-L62)), then the first rate step (`fn optimize_gtr_rate()`, backward passes only, `RateRangeEdge::KeepStart`)
3. For each of the `iterations`: a fit, then a rate step. Each fit counts with the forward messages and the root state of the update of step 2, and with the backward messages of the last rate step, as v0's `infer_gtr()` reuses its existing marginal reconstruction
4. The sampling-bias scale, then a final marginal update under the final model

`fn refine_gtr()` logs the log likelihood at the start and at the end. With `RefinementRate::Keep` (ancestral `--gtr-iterations`), every fit counts with the messages of the start-model update and no rate step runs ([kb/issues/M-ancestral-gtr-iterations-refit-against-frozen-messages.md](../issues/M-ancestral-gtr-iterations-refit-against-frozen-messages.md)).

Pseudo-count smoothing of the initial pi and root state uniform-threshold filtering are opt-in flags (`--smooth-initial-pi`, `--filter-uninformative-root`) that default to v0 behavior. Remaining parity differences are tracked in [kb/issues/H-mugration-iterative-gtr.md](../issues/H-mugration-iterative-gtr.md).

---

## Confidence Profiles

The `profile` of each `StatePosterior` is the full posterior distribution `profile[n_states]` of a node, in the order of the state alphabet. This enables uncertainty quantification for trait assignments and identification of ambiguous nodes (flat profiles, high entropy). The confidence CSV writes the full profile; augur node data and Auspice write the states above 0.001 (`fn build_confidence_map()`) and the entropy.

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

| Test Type          | Location                                                                                                                                                  | Coverage                                                                                                                                  |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Golden master      | [`mugration/__tests__/test_gm_mugration.rs`](../../packages/treetime/src/mugration/__tests__/test_gm_mugration.rs)                                        | States and trait assignments of every node, confidence profiles; nodes named as v0 names them (ladderize, then `NODE_%07d` in preorder) |
| Pipeline           | [`mugration/__tests__/test_run.rs`](../../packages/treetime/src/mugration/__tests__/test_run.rs)                                                          | Weight coverage, `pi`, pseudo-counts, refinement, sampling-bias scale, single-state rejection, leaf without metadata                      |
| Partition          | [`mugration/__tests__/test_discrete_marginal.rs`](../../packages/treetime/src/mugration/__tests__/test_discrete_marginal.rs)                              | Leaf observations and warnings, more than 128 named states, backward and forward passes, full run                                         |
| Transition counts  | [`discrete/__tests__/test_transition_counts.rs`](../../packages/treetime/src/partition/marginal/discrete/__tests__/test_transition_counts.rs)              | Discrete transition counts and node posteriors with 2, 3 and 5 states against a hand-enumerated F81 oracle; root filter                   |
| State posterior    | [`discrete/__tests__/test_state_posterior.rs`](../../packages/treetime/src/partition/marginal/discrete/__tests__/test_state_posterior.rs)                  | Most likely state, ties, entropy                                                                                                          |
| State alphabet     | [`alphabet/__tests__/test_state_alphabet.rs`](../../packages/treetime/src/alphabet/__tests__/test_state_alphabet.rs)                                      | Order, lookup, rejection of fewer than 2 states                                                                                           |
| Refinement         | [`gtr/__tests__/test_refinement.rs`](../../packages/treetime/src/gtr/__tests__/test_refinement.rs)                                                        | Keep and optimize rate policies, messages that each fit counts with, sampling-bias scale                                                  |
| Brent              | inline tests in [`gtr/brent_bracketed.rs`](../../packages/treetime/src/gtr/brent_bracketed.rs)                                                            | Pure math optimizer tests                                                                                                                 |

### Test Gaps

1. Missing data in golden masters: no golden-master dataset has a leaf with an unknown trait, so the v0 parity of imputed leaves is not tested (see "Leaf evidence" in [kb/issues/H-mugration-iterative-gtr.md](../issues/H-mugration-iterative-gtr.md))
2. Confidence parity: four golden-master confidence cases are ignored ([kb/issues/H-mugration-iterative-gtr.md](../issues/H-mugration-iterative-gtr.md))

---

## Known Issues

- [kb/issues/H-mugration-iterative-gtr.md](../issues/H-mugration-iterative-gtr.md): trait assignments match v0 on all seven golden-master datasets; the confidence profiles of four datasets differ from v0 by 3.0e-10 to 2.5e-8
- [kb/issues/M-mugration-gtr-rate-restore-inconsistency.md](../issues/M-mugration-gtr-rate-restore-inconsistency.md): rate restore without a bracket
- [kb/issues/N-mugration-tied-state-follows-rounding.md](../issues/N-mugration-tied-state-follows-rounding.md): a tied most likely state follows rounding

---

## References

- <a id="ref-0"></a>Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042
- <a id="ref-1"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-1)
- <a id="ref-2"></a>Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. [↩](#cite-2)
- <a id="ref-3"></a>Ishikawa, Sohta A., Anna Zhukova, Wataru Iwasaki, and Olivier Gascuel. 2019. "A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios." _Molecular Biology and Evolution_ 36(9):2069-2085. https://doi.org/10.1093/molbev/msz131 [↩](#cite-3)
- <a id="ref-4"></a>Vaiente, Mara A., and Matthew Scotch. 2020. "Going Back to the Roots: Evaluating Bayesian Phylogeographic Models with Discrete Trait Uncertainty." _Infection, Genetics and Evolution_ 85:104501. https://doi.org/10.1016/j.meegid.2020.104501 [↩](#cite-4)
- <a id="ref-5"></a>Lemey, Philippe, Andrew Rambaut, Alexei J. Drummond, and Marc A. Suchard. 2009. "Bayesian Phylogeography Finds Its Roots." _PLoS Computational Biology_ 5(9):e1000520. https://doi.org/10.1371/journal.pcbi.1000520
- <a id="ref-6"></a>Edwards, Ceiridwen J., Marc A. Suchard, Philippe Lemey, et al. 2011. "Ancient Hybridization and an Irish Origin for the Modern Polar Bear Matriline." _Current Biology_ 21(15):1251-1258. https://doi.org/10.1016/j.cub.2011.05.058
