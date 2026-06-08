# Standalone branch length optimization command

v1 exposes branch length optimization as a standalone CLI command `treetime optimize`. v0 performs branch length optimization only as an internal step within `optimize_tree()` during ancestral reconstruction and timetree inference, with no user-facing entry point.

## What v0 does

v0 optimizes branch lengths inside `TreeAnc.optimize_tree()` ([packages/legacy/treetime/treetime/treeanc.py#L1384-L1473](../../packages/legacy/treetime/treetime/treeanc.py#L1384-L1473)), which alternates between ancestral reconstruction and per-branch optimization. Users interact with it indirectly through the `ancestral` and default timetree commands. The branch optimization method, convergence parameters, and damping factor are hardcoded or controlled by internal flags not exposed on the CLI.

The per-branch optimizer `GTR.optimal_t_compressed()` ([packages/legacy/treetime/treetime/gtr.py#L816-L907](../../packages/legacy/treetime/treetime/gtr.py#L816-L907)) uses Brent's method via `scipy.optimize.minimize_scalar`. The outer loop applies exponential damping (factor 0.75) to branch length updates for stability.

There is no way for a v0 user to optimize branch lengths on a tree without also performing ancestral reconstruction or timetree inference.

## What v1 does

v1 adds a dedicated `treetime optimize` subcommand ([packages/treetime/src/commands/optimize/run.rs#L34-L145](../../packages/treetime/src/commands/optimize/run.rs#L34-L145)) that optimizes branch lengths given a tree and alignment. The command accepts user-facing parameters for convergence control.

### CLI arguments

| Argument                   | Default                                  | Purpose                                                    |
| -------------------------- | ---------------------------------------- | ---------------------------------------------------------- |
| `--tree`                   | (required)                               | Input tree in Newick/Nexus/Phylip                          |
| `--aln`                    | (required)                               | FASTA alignment (supports compression)                     |
| `--model`                  | `infer`                                  | GTR model (JC69, K80, F81, HKY85, T92, TN93, JTT92, infer) |
| `--alphabet`               | auto                                     | Nucleotide or amino acid alphabet                          |
| `--dense`                  | auto                                     | Force dense sequence representation                        |
| `--max-iter`               | 10                                       | Maximum outer iterations                                   |
| `--dp`                     | 0.1                                      | Convergence tolerance (log-likelihood delta)               |
| `--outdir`                 | (required)                               | Output directory                                           |
| `--output-augur-node-data` | `<outdir>/optimize.augur-node-data.json` | Path for the augur-compatible node data JSON               |

### Algorithm

The command runs an iterative cycle of marginal ancestral reconstruction and per-branch optimization:

1. Read alignment and tree. Create both sparse and dense partitions from the same alignment.
2. Run Fitch compression on sparse partitions. Initialize marginal reconstruction (backward and forward passes) on both partition types.
3. Compute initial branch lengths as Hamming distance divided by sequence length across all partitions.
4. Iterate until convergence or `max-iter`:
   - Run marginal reconstruction on both sparse and dense partitions, computing total log-likelihood.
   - If `|log_lh - log_lh_prev| < dp`, stop.
   - Optimize each branch length independently using Newton-Raphson with grid search fallback (documented separately in [optimize-newton-raphson-per-edge.md](optimize-newton-raphson-per-edge.md)).
5. Write `annotated_tree.nwk`, `annotated_tree.nexus`, `gtr.json`, and `optimize.augur-node-data.json` to the output directory.

### Augur node data JSON output

The command writes an augur-compatible node data JSON ([packages/treetime/src/commands/optimize/augur_node_data.rs](../../packages/treetime/src/commands/optimize/augur_node_data.rs)), consumed by `augur export v2 --node-data`. There is no augur "optimize" command, so the output contract is `augur refine` run WITHOUT `--timetree`: top-level `generated_by`, `alignment`, and `input_tree`, plus per-node `branch_length`. Augur's non-timetree path sets `attributes = ['branch_length', 'confidence']` and skips the clock, date, and `num_date_confidence` fields ([augur/refine.py](https://github.com/nextstrain/augur/blob/024292af6daf/augur/refine.py)).

The one content difference from augur's non-timetree refine: augur reads `branch_length` from the unmodified input tree (it only instantiates `TreeAnc` to name internal nodes), whereas `optimize` writes the ML-optimized branch length into the same field. The file shape is identical; the values are v1's optimized divergences (substitutions per site). `augur export v2`'s `node_div()` consumes `branch_length` for cumulative divergence when `mutation_length` is absent.

The types are reused from the `util-augur-node-data-json` crate (`AugurNodeDataJsonRefine`), shared with the timetree command. `confidence` is omitted because v1's Newick reader does not parse input-tree branch support values ([kb/issues/N-timetree-node-data-confidence-not-emitted.md](../issues/N-timetree-node-data-confidence-not-emitted.md)). Both commands support `--divergence-units=mutations` to emit integer mutation counts instead of float subs/site values.

The per-branch optimization uses the GTR eigensystem decomposition to precompute branch-length-independent coefficients `k_c = (msg_child . V)_c * (msg_parent . V_inv^T)_c` for each alignment position. These coefficients allow computing log-likelihood, first derivative, and second derivative from the same cached exponentials. Dense partitions contribute per-position coefficients; sparse partitions group invariant positions by state pair with multiplicity counts ([packages/treetime/src/partition/optimize_sparse.rs#L36-L119](../../packages/treetime/src/partition/optimize_sparse.rs#L36-L119)).

The same optimization machinery is reused by the timetree command to compute branch length likelihood distributions on a grid ([packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L31-L63](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs#L31-L63)).

## Scientific background

### Maximum likelihood branch length estimation

Branch lengths in a phylogenetic tree represent expected substitutions per site along each edge. Given a fixed topology and a sequence alignment, maximum likelihood (ML) estimation finds the branch lengths that maximize the probability of observing the alignment data under a substitution model.

Felsenstein's pruning algorithm (Felsenstein, 1981) computes the likelihood of observing an alignment given a tree by recursing from leaves to root. At each internal node k with children i, j:

```
w_k(X) = [sum_Y P(X->Y, t_i) * w_i(Y)] * [sum_Z P(X->Z, t_j) * w_j(Z)]
```

where `P(X->Y, t)` is the transition probability from state X to Y over branch length t, computed from the substitution model rate matrix Q as `P(t) = exp(Qt)`.

The total log-likelihood is the sum over all alignment sites of the log of the root partial likelihood weighted by equilibrium frequencies.

### Eigendecomposition for efficient optimization

The GTR rate matrix Q can be decomposed as `Q = V * diag(lambda) * V^{-1}` where `lambda` are eigenvalues and V is the eigenvector matrix. This transforms the matrix exponential into:

```
P(t) = V * diag(exp(lambda_c * t)) * V^{-1}
```

V and V^{-1} are computed once per model. For each branch length t, only the diagonal `exp(lambda_c * t)` changes. The derivatives follow directly:

```
dP/dt  = V * diag(lambda_c * exp(lambda_c * t)) * V^{-1}
d^2P/dt^2 = V * diag(lambda_c^2 * exp(lambda_c * t)) * V^{-1}
```

This makes Newton-Raphson practical: likelihood, gradient, and curvature share the same exponential evaluations.

### Coordinate-wise optimization

Standard practice in phylogenetic ML software (RAxML, IQ-TREE, PhyML) is coordinate descent: optimize one branch length at a time, cycling through all branches, repeating until convergence. Clancy, Lyu, and Roch (arXiv 2507.22038) proved that for the two-state symmetric model on balanced binary trees, coordinate-wise ML converges exponentially fast to the MLE under standard regularity conditions. Each branch optimization is a one-dimensional problem with at most one stationary point under JC69 and F81 (Dinh and Matsen, arXiv 1507.03647).

## Why v1 adds this command

**Separation of concerns.** v0 bundles branch length optimization with ancestral reconstruction and timetree inference. Exposing optimization as a standalone command enables workflows where a user wants to refine branch lengths on an existing tree without performing temporal inference or ancestral state assignment.

**Alignment quality assessment.** Comparing the log-likelihood before and after optimization reveals how well the substitution model fits the data. A large likelihood improvement indicates the input tree's branch lengths are inconsistent with the alignment under the chosen model.

**Pipeline composability.** The command reads standard tree and alignment files and writes annotated trees. It slots into pipelines where a tree from one tool (parsimony, neighbor-joining, external ML inference) needs branch length refinement under a specific substitution model before downstream analysis.

## Differences from v0's internal optimization

| Aspect              | v0 (internal)                            | v1 (standalone command)                       |
| ------------------- | ---------------------------------------- | --------------------------------------------- |
| Exposure            | Internal method within `optimize_tree()` | Standalone CLI subcommand                     |
| Convergence control | Hardcoded parameters                     | User-specified `--max-iter` and `--dp`        |
| Outer loop damping  | Exponential decay (factor 0.75)          | Exponential decay (`--damping`, default 0.75) |
| Per-branch method   | Brent's method (derivative-free)         | Newton-Raphson + grid fallback                |
| Partition support   | Single representation                    | Mixed dense + sparse partitions               |
| Model selection     | GTR inferred internally                  | User-specified via `--model`                  |

The per-branch optimization method difference (Newton-Raphson vs Brent) is documented separately in [optimize-newton-raphson-per-edge.md](optimize-newton-raphson-per-edge.md).

## Practical impact

Users can run `treetime optimize --tree=tree.nwk --aln=aln.fasta --outdir=out/` to refine branch lengths without performing ancestral reconstruction or timetree inference. The output tree retains the input topology with updated branch lengths.

v1 applies the same outer-loop exponential damping as v0 (`--damping`, default 0.75) to prevent oscillation between iterations for datasets where the likelihood surface is poorly conditioned.

## References

- Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359
- Dinh, Vu, and Frederick A. Matsen IV. 2015. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." arXiv:1507.03647. https://arxiv.org/abs/1507.03647
- Clancy, Killian M., Hanlin Lyu, and Sebastien Roch. 2025. "Sample Complexity of Branch-Length Estimation by Maximum Likelihood." arXiv:2507.22038. https://arxiv.org/abs/2507.22038
- Stamatakis, Alexandros. 2006. "RAxML-VI-HPC: Maximum Likelihood-Based Phylogenetic Analyses with Thousands of Taxa and Mixed Models." _Bioinformatics_ 22(21):2688-2690. https://doi.org/10.1093/bioinformatics/btl446
- Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Molecular Biology and Evolution_ 32(1):268-274. https://doi.org/10.1093/molbev/msu300
- Guindon, Stephane, and Olivier Gascuel. 2003. "A Simple, Fast, and Accurate Algorithm to Estimate Large Phylogenies by Maximum Likelihood." _Systematic Biology_ 52(5):696-704. https://doi.org/10.1080/10635150390235520
