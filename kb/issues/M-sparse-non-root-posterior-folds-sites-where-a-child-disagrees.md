# Sparse posterior of a non-root node folds sites where a child disagrees

Below the root, the sparse posterior of a node differs from the dense posterior at positions where a child disagrees with the node: the child subtree has no data there, or a substitution on the child edge gives the child another state. The forward pass folds such a position into the shared profile of its state once the posterior exceeds $1 - 10^{-4}$, and the shared profile assumes that every child carries the node state. The root keeps the same positions explicit and agrees with dense.

## Mechanism

The backward pass keeps a position explicit in a node's message to its parent when a child disagrees. The candidate site table of `fn combine_child_messages()` ([packages/treetime/src/partition/marginal/sparse/backward.rs](../../packages/treetime/src/partition/marginal/sparse/backward.rs)) gives the child a gap or unknown reference at its no-data sites and the query state of a substitution on its edge, so `fn combine_messages()` ([packages/treetime/src/partition/marginal/sparse/message.rs](../../packages/treetime/src/partition/marginal/sparse/message.rs)) sees a message that differs from the node state, and the agreement condition of the fold fails. The position then folds only when its profile equals the shared profile of its state bit for bit (`fn is_shared_profile()`). At the root, this product is the posterior.

Below the root, `fn SparseForwardPass::edge()` ([packages/treetime/src/partition/marginal/sparse/forward.rs](../../packages/treetime/src/partition/marginal/sparse/forward.rs)) rebuilds the posterior from two messages: the message from the parent and the node's own message to its parent. The site table of this product has no references of the children. The explicit entry of the node's message carries the node state, and the message from the parent agrees with it, so the agreement condition holds and the position folds when its posterior exceeds $1 - 10^{-4}$. `fn compute_msg_to_child()` then divides the shared profile, instead of the posterior, by the message of each child, so the outside messages of the children and their branch-length objectives lose the difference.

## History of the threshold

- `5c59ff2df` (port of the sparse prototype): a position stays explicit only when its subleading states exceed $\epsilon = 10^{-6}$; all others fold into the shared profile. The fold keeps the number of explicit sites small
- `b7a0f7ba9` sets $\epsilon$ to $10^{-4}$, after $10^{-3}$
- `ff26255dd` (sparse branch-length optimization) adds the agreement condition: a position where a message disagrees with the node state never folds by the threshold
- `303086b14` adds the exact fold: a position whose profile equals the shared profile of its state folds also where a message disagrees

The forward pass has no references of the children, so the agreement condition of `ff26255dd` cannot see a disagreeing child there. The threshold guards no property that an exact posterior breaks.

## Evidence

The dense representation is the oracle. `test_root_polytomy_sparse_matches_dense` in [packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs) compares the root log-likelihood, the internal-node posteriors, the change of every edge objective at three branch lengths, and the optimal branch length of every edge with dense, for a polytomy at the root and one level below it (node `P` at branch length 0.02, next to a tip). Each tip differs from a random root sequence by one substitution. Largest differences, sparse minus dense:

| Case | Tips, branch length, `N` runs per tip | Edge objective change | Optimal branch length | Posterior | Root log-likelihood |
| --- | --- | ---: | ---: | ---: | ---: |
| `long_branches_at_root` | 8, 0.8, 1 | 4.3e-14 | -1.3e-13 | 2.2e-16 | 0 |
| `long_branches_below_root` | 8, 0.8, 1 | -0.043 (edge `c6`, t = 0.001) | 4.4e-5 (edge `c6`) | 7.3e-5 (`P`, position 21: dense 0.999922, sparse 0.999995) | 2.3e-13 |
| `long_branches_with_data_at_root` | 12, 0.8, 0 | 1.1e-13 | -2.0e-12 | -4.4e-16 | -1.1e-13 |
| `long_branches_with_data_below_root` | 12, 0.8, 0 | -1.5e-3 (edge `c8`, t = 0.001) | 1.5e-6 (edge `c8`) | 1.9e-6 (`P`, position 52: dense 0.9999981, sparse 0.99999995) | -4.5e-13 |
| `short_branches_at_root` | 12, 0.01, 2 | 7.1e-15 | 2.9e-14 | 2.2e-16 | 2.8e-14 |
| `short_branches_below_root` | 12, 0.01, 2 | 2.8e-14 | 5.0e-14 | 2.2e-16 | 1.1e-13 |

The case without `N` runs shows that substitution sites fold in the same way as sites without data. The error grows with the branch lengths of the children: at branch length 0.01 it is below round-off. With the exact rule of the fix direction, all six cases agree with dense within 2e-12.

On the public dataset `data/flu/h3n2/500`, `optimize` with a fixed topology gives these branch lengths, sparse against dense over 917 edges:

| Sparse run | Largest difference | Edges that differ by more than 1e-7 |
| --- | ---: | ---: |
| `rust` | 3.4e-4 (tip `A/Morioka/22/2006`: sparse 3.68e-4, dense 7.10e-4) | 398 |
| Exact rule | 1.0e-6 | 23 |

On `data/mpox/clade-ii/500`, which has short branches, the exact rule moves the `optimize` branch lengths by at most 1.6e-9.

## Effect on a production build

On a production mpox B.1 build (4,946 tips, 197,209 sites, a root polytomy of 657 children, `timetree` with marginal branch lengths), the exact rule changes:

- `log_lh_seq` by +26.6 nats (-341165.61 against -341139.05 in the last iteration)
- 1,147 node dates by more than 1e-9 years: 99th percentile 0.0012 years, maximum 0.008 years

In one forward pass, the fold removes about 5.0M explicit posterior entries, summed over the nodes. The exact rule keeps them; about 4.4k of them differ from the shared profile by 1e-9 or more (2 by more than 0.1, about 2.4k by more than 1e-6), about 1.0M by 1e-12 to 1e-9, 0.95M by 1e-16 to 1e-12, and 3.06M by less than 1e-16.

## Reproduction

The below-root cases are ignored. Run them on `rust`; the two below-root cases with long branches fail:

```bash
./dev/docker/run just t -E 'package(treetime) & test(test_root_polytomy)' --run-ignored all
```

On a public dataset, run `optimize` with a fixed topology in both representations and compare the `branch_length` of each node in `optimize.node-data.json`:

```bash
./dev/docker/run just r treetime optimize --representation=sparse --tree=data/flu/h3n2/500/tree.nwk --alignment=data/flu/h3n2/500/aln.fasta.xz --no-collapse-short-branches --no-merge-shared-mutations --no-hoist-reversions --output-all=tmp/fold/sparse
./dev/docker/run just r treetime optimize --representation=dense --tree=data/flu/h3n2/500/tree.nwk --alignment=data/flu/h3n2/500/aln.fasta.xz --no-collapse-short-branches --output-all=tmp/fold/dense
```

## Fix direction

The exact rule: the node's message to its parent enters the forward product under a fold rule that lets its explicit entries fold only when they equal the shared profile bit for bit, while the message from the parent keeps the threshold fold. The backward pass already applied the fold rule at the node with the references of all children, so the posterior of a node below the root then keeps the positions that the root keeps. Commit `7b4e0ba78` on branch `worktree/fix/sparse-posterior-below-root` implements it (not landed) and enables the ignored cases.

The exact rule doubles the run time of the production mpox B.1 build: 552 s against 270 s, one thread, both runs under the same load. In one forward pass, the explicit posterior entries grow from 15.3k to 5.02M and the outside-message entries from 5.05M to 33.0M, although only 62k and 636k of them are bit-distinct. [kb/proposals/sparse-shared-profiles-for-sites-without-data.md](../proposals/sparse-shared-profiles-for-sites-without-data.md) describes a sharing mechanism that removes these copies without a tolerance.

> [!IMPORTANT]
> **Decision required.** `rust` keeps the threshold fold below the root until the decision on the proposal. Options: take the exact rule with single-site entries (sparse matches dense, the production build takes twice as long), take it together with the class-keyed shared profiles of the proposal (exact, not implemented), or keep the fold.

## Related issues

- [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md): names the fold of resolved positions as one lead for the internal-node residues that differ between the representations
- [M-sparse-dense-likelihood-residual-without-missing-data.md](M-sparse-dense-likelihood-residual-without-missing-data.md): a likelihood difference without missing data, where substitution sites fold below the root
