# Initial branch length estimation: algorithms, context, and design space

## The optimization problem

Branch length optimization finds the evolutionary distance `t` for each edge that maximizes the probability of observing the data. For a single edge, the per-site likelihood under a GTR substitution model is:

```
L_i(t) = sum_{a,b} pp_i(a) * exp(Q*t)_{ab} * pc_i(b)
```

where `pp_i` and `pc_i` are the state distributions at position `i` for the parent and child endpoints, `Q` is the rate matrix, and `exp(Qt)` is the transition probability matrix.

The total log-likelihood sums over all positions:

```
log L(t) = sum_i log L_i(t)
```

Finding the ML branch length requires a scalar optimizer (Brent, Newton-Raphson, grid search) applied to each edge independently, with all other branch lengths and ancestral state distributions held fixed. The optimizer needs a starting point.

## What the starting point is for

The starting point seeds the optimizer. It does not need to be accurate - it needs to be in the basin of attraction of the correct optimum.

**Brent's method** (v0) takes a bracket `[a, m, b]` where `m` is the midpoint guess. Brent searches the interval by combining golden section with parabolic interpolation. The midpoint determines where Brent starts looking, influencing convergence speed but not the final result (for unimodal surfaces).

**Newton-Raphson** (v1) starts from the current branch length and iterates `t_{n+1} = t_n - L'(t_n)/L''(t_n)`. For concave surfaces, convergence is quadratic from any starting point. The starting point determines iteration count and whether the method overshoots into a non-concave region.

**Grid search** (v1 fallback) evaluates likelihood at 100 equally spaced points on `[0.1/L, 1.5*t + 1/L]` where `t` is the current branch length and `L` is the sequence length. The search range depends on `t`, so an underestimated starting point narrows the range.

## The Hamming distance family

All approaches estimate the initial branch length as `d/n` where `d` counts (or estimates) the number of differing positions and `n` is the effective alignment length. They differ in how `d` is computed.

### Hard Hamming on discrete sequences

```
d = count(parent_state[i] != child_state[i])  for non-gap i
```

Each position contributes exactly 0 or 1. Parent and child states come from a discrete reconstruction (parsimony, joint ML, or MAP of marginal posterior).

This is v0's `profiles=False` path. The `state_pair` method at [packages/legacy/treetime/treetime/gtr.py#L631](../../packages/legacy/treetime/treetime/gtr.py#L631) compresses the alignment into unique (parent, child) state pairs with multiplicities, then:

```python
hamming_distance = np.sum(multiplicity[seq_pair[:, 1] != seq_pair[:, 0]]) / np.sum(multiplicity)
```

### Soft Hamming (profile overlap)

```
d = sum_i (1 - dot(pp_i, pc_i))    for non-gap i
```

where `pp_i` and `pc_i` are probability vectors over states at position `i`. Each position contributes a value in `[0, 1]`.

The dot product `dot(pp_i, pc_i) = sum_a pp_i(a) * pc_i(a)` is the Bhattacharyya coefficient: the probability that independent draws from the two distributions yield the same state. The complement `1 - dot(pp_i, pc_i)` measures the probability of disagreement.

Properties:

- When profiles are one-hot (sharp): `dot = 1` if same state, `dot = 0` if different. Reduces to hard Hamming.
- When both uniform over `n` states: `dot = n * (1/n)^2 = 1/n`. Contribution = `1 - 1/n` = `0.75` for nucleotides.
- Continuous and differentiable with respect to profile changes.
- Monotonically related to the Hellinger distance: `H^2 = 1 - sum sqrt(pp * pc)`.

This is v0's `profiles=True` path:

```python
hamming_distance = 1 - np.sum(multiplicity * np.sum(seq_pair[0] * seq_pair[1], axis=1)) / np.sum(multiplicity)
```

### The critical question: profiles of what?

The soft Hamming formula is well-defined for any pair of probability distributions, but its meaning depends on what the distributions represent. Three choices have been considered in this project:

**Choice B: Both distributions at the child node** (v0)

v0's `marginal_branch_profile(node)` at [packages/legacy/treetime/treetime/treeanc.py#L1122-L1146](../../packages/legacy/treetime/treetime/treeanc.py#L1122-L1146) returns:

- `pp = node.marginal_outgroup_LH` - outgroup evidence, evaluated at the child position
- `pc = node.marginal_subtree_LH` - subtree evidence, evaluated at the child position

Both distributions describe beliefs about the **same node** (the child), derived from two independent sources of evidence. The outgroup message `pp` was propagated through the transition matrix from the parent side to the child position during the forward pass of marginal inference.

`dot(pp_i, pc_i)` = probability that the two independent evidence sources agree on the child's state at position `i`. When they disagree strongly, it indicates the branch carries a substitution at that position. This has a clean probabilistic interpretation.

**Choice C: Messages at opposite ends of the edge** (removed v1 override)

The removed code at `marginal_dense.rs:131-152` used:

- `msg_to_parent.dis` - subtree evidence at the **child** end
- `msg_to_child.dis` - outgroup evidence at the **parent** end

These distributions describe beliefs about **different nodes**. `dot(msg_to_parent_i, msg_to_child_i)` does not equal the probability that parent and child share a state, because the transition matrix connecting the two positions is not accounted for.

Example: `msg_to_parent = [0.9, 0.1, 0, 0]` (child is probably A), `msg_to_child = [0.6, 0.4, 0, 0]` (parent is probably A). The dot product is `0.54 + 0.04 = 0.58`. But the actual probability that parent and child share a state depends on the branch length and transition matrix, which this dot product ignores.

This formula produces continuous values in the right range but is not grounded in a valid probabilistic measure.

**Choice E: Posterior distributions at parent and child nodes** (proposed, not implemented)

Use the full marginal posteriors at the parent and child endpoints:

- `parent_posterior = nodes[parent_key].profile.dis` - full posterior at parent
- `child_posterior = nodes[child_key].profile.dis` - full posterior at child

Each posterior incorporates all tree evidence (subtree + outgroup) at its respective node. `dot(parent_posterior_i, child_posterior_i)` = probability that independent draws from the parent and child posteriors yield the same state at position `i`.

This is a valid probabilistic measure, but note: the two distributions describe different nodes (parent and child), not the same node from two evidence sources (like v0). The interpretation is: "how similar are our best estimates of the parent and child states?" rather than "how much do the two evidence sources disagree about the child?"

## Comparison table

| Approach                         | Formula                         | Data source                | Position               | Continuity       | Probabilistic interpretation       |
| -------------------------------- | ------------------------------- | -------------------------- | ---------------------- | ---------------- | ---------------------------------- |
| A. Hard Hamming (discrete)       | `I(s_p != s_c)`                 | Discrete sequences         | parent + child         | Discrete {0,1}   | Indicator of state mismatch        |
| B. Soft Hamming (child-end)      | `1 - dot(pp, pc)`               | Marginal profiles at child | both at child          | Continuous [0,1] | P(outgroup disagrees with subtree) |
| C. Soft Hamming (edge messages)  | `1 - dot(msg_p, msg_c)`         | Per-edge partial messages  | parent-end + child-end | Continuous [0,1] | None (cross-position overlap)      |
| D. Hard Hamming (MAP posteriors) | `I(argmax(π_p) != argmax(π_c))` | Full node posteriors       | parent + child         | Discrete {0,1}   | Indicator of MAP state mismatch    |
| E. Soft Hamming (posteriors)     | `1 - dot(π_p, π_c)`             | Full node posteriors       | parent + child         | Continuous [0,1] | P(posterior draws disagree)        |

### Numerical examples (4-state nucleotide alphabet)

```
Position type                    A      B      C      D      E
─────────────────────────────────────────────────────────────────
Sharp, same state               0.00   0.00   0.00   0.00   0.00
Sharp, different states         1.00   1.00   1.00   1.00   1.00
Both uniform                    0.00   0.75   0.75   0.00   0.75
Weak same ([.7,.1,.1,.1] both)  0.00   0.48   varies 0.00   0.48
Weak diff ([.7,.1,.1,.1] vs     1.00   0.84   varies 1.00   0.84
  [.1,.7,.1,.1])
```

The critical divergence is at uncertain positions. A and D give 0 (argmax agrees), while B and E give ~0.75 (profiles have low overlap). C gives model-dependent values with no clean relationship to the others.

## How v0 uses the Hamming distance

v0's `optimize_tree_marginal()` at [packages/legacy/treetime/treetime/treeanc.py#L1297-L1360](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360) loops:

```
for each iteration:
    for each non-root node:
        pp, pc = marginal_branch_profile(node)                    # child-end profiles
        new_bl = optimal_t_compressed((pp,pc), mult, profiles=True) # Brent optimization
        node.branch_length = damped(new_bl, old_bl, iteration)    # exponential damping
    infer_ancestral_sequences(marginal=True)                      # recompute profiles
```

Inside `optimal_t_compressed(profiles=True)` at [packages/legacy/treetime/treetime/gtr.py#L816-L920](../../packages/legacy/treetime/treetime/gtr.py#L816-L920):

1. Compute soft Hamming: `h = 1 - sum(mult * sum(pp * pc)) / sum(mult)`
2. Use as Brent bracket midpoint: `bracket = [-sqrt(MAX), sqrt(h), sqrt(MAX)]`
3. Brent optimizes `prob_t_profiles` (full profile likelihood) in sqrt(t) space
4. If Brent fails: fall back to `h` as the branch length

The profile likelihood `prob_t_profiles` at [packages/legacy/treetime/treetime/gtr.py#L922-L963](../../packages/legacy/treetime/treetime/gtr.py#L922-L963) is:

```python
res = np.einsum('ai,ij,aj->a', pc, Qt, pp)   # per-site likelihood
logP = sum(mult * log(res))                    # total log-likelihood
```

This is the exact same eigendecomposition-based likelihood that v1 uses (different notation, same math).

Key observations:

- v0's soft Hamming and Brent likelihood use the **same data** (`pp, pc` profiles)
- The Hamming distance is a rough approximation to the Brent-optimal branch length
- v0 also adds a regularization penalty `exp(t^4/10000)` to prevent runaway for profile mode
- Damping (`0.75` default) blends old and new branch lengths to stabilize the outer loop

## How v1 uses the initial guess

v1's `run_optimize()` at [packages/treetime/src/commands/optimize/run.rs#L34-L150](../../packages/treetime/src/commands/optimize/run.rs#L34-L150):

```
update_marginal(sparse)           # sparse marginal inference
initialize_marginal(dense)         # dense initialization
update_marginal(dense)             # dense marginal (dummy JC69)
replace GTR with real model
update_marginal(dense)             # dense marginal (real GTR)
initial_guess_mixed(partitions)    # <--- sets initial branch lengths
for each iteration:
    update_marginal(sparse + dense)
    run_optimize_mixed(partitions)  # Newton-Raphson per edge
```

`initial_guess_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L293-L325](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L293-L325):

```rust
for each edge:
    sub_count = sum over partitions: edge_subs().len()
    effective_length = sum over partitions: edge_effective_length()
    branch_length = sub_count / effective_length
```

`run_optimize_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L216-L289](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L216-L289):

```
for each edge:
    coefficients = msg_to_child.dot(V) * msg_to_parent.dot(V_inv.T)   # eigenvalue-space
    if is_zero_branch_optimal: set 0, continue
    Newton-Raphson using L, L', L'' from coefficients
    Fallback: grid search on [0.1/L, 1.5*t + 1/L]
```

Key observations:

- v1's initial guess and optimizer use **different data**: initial guess reads node posteriors via `edge_subs()`, optimizer reads per-edge messages
- The initial guess runs once before the loop; the optimizer runs every iteration
- v1 has no damping in the outer loop (known issue: `M-optimize-oscillation-no-damping.md`)
- v1 has no regularization penalty (known issue - not explicitly tracked)

## Design space analysis

### For the initial guess specifically

The initial guess is consumed by Newton-Raphson on the first optimization iteration. Its quality matters for:

1. **Convergence speed**: closer starting point = fewer Newton iterations
2. **Convergence reliability**: for multimodal surfaces (K80+ models per Dinh & Matsen 2015), the starting point determines which optimum is found
3. **Grid search range**: when Newton fails, the fallback grid is `[0.1/L, 1.5*t + 1/L]`. If `t=0` (hard Hamming at a fully uncertain edge), the range narrows to `[0.1/L, 1/L]` - covering only zero-to-one-mutation branch lengths

### Option D (current: hard Hamming on MAP posteriors)

Advantages:

- Simple: no dot products, just argmax comparison
- Uses correct data source (full posteriors)
- Consistent with sparse behavior
- Explicitly requested by the user

Disadvantages:

- Loses uncertainty information: two positions with identical argmax but very different posterior shapes both contribute 0
- Systematically underestimates branch lengths at uncertain edges
- Narrows grid search range for uncertain edges

### Option E (proposed: soft Hamming on node posteriors)

Advantages:

- Continuous: smooth with respect to posterior changes
- Valid probabilistic interpretation: `dot(π_p, π_c)` is the probability of state agreement
- Preserves uncertainty information
- Closer to v0's scientific intent (profile-based estimation)

Disadvantages:

- Different from v0's exact formula (v0 evaluates both profiles at the child, E evaluates at parent and child)
- More expensive: dot product per position per edge
- For very uncertain positions (uniform posteriors), gives 0.75 per position, which can overestimate branch lengths

### Option B (v0-faithful: soft Hamming at child end)

To reproduce v0 exactly, you'd need:

1. The outgroup message at the child node (propagated through expQt from parent)
2. The subtree message at the child node

In v1's architecture, `msg_to_child` is the outgroup message at the **parent** end. To evaluate it at the child end, you'd need to propagate: `pp_at_child = msg_to_child.dot(expQt)`. But computing `expQt` requires the current branch length - which is what we're trying to estimate.

v0 avoids this chicken-and-egg because the forward pass already propagated the outgroup message to the child node (stored as `node.marginal_outgroup_LH`). v1's forward pass stores the propagated outgroup message in the node posterior (`nodes[key].profile.dis`), which combines subtree and outgroup evidence. The separate outgroup-at-child profile is not stored.

Options to implement B:

- Store the propagated outgroup-at-child profile separately during the forward pass (architecture change)
- Recompute it on the fly using the current branch length (chicken-and-egg, but could use tree-file branch length as approximation)
- Accept that E is a reasonable approximation to B for well-resolved trees

### Option "None" (use tree-file branch lengths)

Skip `initial_guess_mixed()` entirely and let Newton start from whatever the tree file provides. This works when:

- The tree was estimated by another tool (RAxML, IQ-TREE) with reasonable branch lengths
- The tree was previously optimized by treetime

This is what v1's timetree command currently does (known issue: `M-timetree-missing-initial-branch-optimization.md`).

## Summary of the current state

After this branch, v1 uses option D (hard Hamming on MAP posteriors). This is a deliberate v0 divergence documented in the Newton-Raphson intentional change doc. The divergence is narrow:

- For well-resolved trees (sharp posteriors): D gives the same result as B, C, or E
- For trees with many uncertain internal edges: D underestimates initial branch lengths compared to B or E
- Newton-Raphson corrects the underestimate during optimization, except when it falls back to the grid search

If future testing reveals convergence problems at uncertain edges (grid search missing optima due to narrow range from D's zero estimates), option E would be the natural fix: replace `edge_subs().len()` with a dot-product loop over node posteriors. This requires no architecture changes and adds ~20 lines to `initial_guess_mixed()`.
