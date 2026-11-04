# Sparse posterior of a non-root node folds sites where a child disagrees

Below the root, the sparse posterior of a node differs from the dense posterior at positions where a child disagrees with the node: the child subtree has no data there, or a substitution on the child edge gives the child another state. The forward pass folds such a position into the <a id="gloss-use-1"></a>shared profile <sup>[1](#gloss-1)</sup> of its state once the posterior exceeds $1 - 10^{-4}$, and the shared profile assumes that every child carries the node state. The root keeps the same positions explicit and agrees with dense.

## Mechanism

The backward pass keeps a position explicit in a node's message to its parent when a child disagrees. The candidate site table of `fn combine_child_messages()` ([packages/treetime/src/partition/marginal/sparse/backward.rs](../../packages/treetime/src/partition/marginal/sparse/backward.rs)) gives the child a gap or unknown reference at its no-data sites and the query state of a substitution on its edge, so `fn combine_messages()` ([packages/treetime/src/partition/marginal/sparse/message.rs](../../packages/treetime/src/partition/marginal/sparse/message.rs)) sees a message that differs from the node state, and the agreement condition of the fold fails. The position then folds only when its profile equals the shared profile of its state bit for bit (`fn is_shared_profile()`). At the root, this product is the posterior.

Below the root, `fn SparseForwardPass::edge()` ([packages/treetime/src/partition/marginal/sparse/forward.rs](../../packages/treetime/src/partition/marginal/sparse/forward.rs)) rebuilds the posterior from two messages: the message from the parent and the node's own message to its parent. The site table of this product has no references of the children. The explicit entry of the node's message carries the node state, and the message from the parent agrees with it, so the agreement condition holds and the position folds when its posterior exceeds $1 - 10^{-4}$. `fn compute_msg_to_child()` then divides the shared profile, instead of the posterior, by the message of each child, so the <a id="gloss-use-2"></a>outside messages <sup>[2](#gloss-2)</sup> of the children and their branch-length objectives lose the difference.

## History of the threshold

- `5c59ff2df` (port of the sparse prototype): a position stays explicit only when its subleading states exceed $\epsilon = 10^{-6}$; all others fold into the shared profile. The fold keeps the number of explicit sites small
- `b7a0f7ba9` sets $\epsilon$ to $10^{-4}$, after $10^{-3}$
- `ff26255dd` (sparse branch-length optimization) adds the agreement condition: a position where a message disagrees with the node state never folds by the threshold
- `303086b14` adds the exact fold: a position whose profile equals the shared profile of its state folds also where a message disagrees

The forward pass has no references of the children, so the agreement condition of `ff26255dd` cannot see a disagreeing child there. The threshold guards no property that an exact posterior breaks.

## Evidence

The dense representation is the oracle. `test_root_polytomy_sparse_matches_dense` in [packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs) compares the root log-likelihood, the internal-node posteriors, the change of every <a id="gloss-use-3"></a>edge objective <sup>[3](#gloss-3)</sup> at three branch lengths, and the optimal branch length of every edge with dense, for a polytomy at the root and one level below it (node `P` at branch length 0.02, next to a tip). Each tip differs from a random root sequence by one substitution. Largest differences, sparse minus dense:

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

One `timetree` run with the threshold fold of `rust` against one run with the exact rule of the fix direction, same inputs and same host, on a production mpox B.1 build: 4,946 tips, 6,396 nodes after the collapse of short branches, 197,209 sites, a root polytomy of 657 children, marginal branch lengths, and a fixed clock rate of $5.7 \cdot 10^{-5}$ substitutions per site per year. Dense mode was not run on this build; the exact rule matches dense within 2e-12 on the test cases above. Log-likelihoods of the last iteration, in <a id="gloss-use-4"></a>nats <sup>[4](#gloss-4)</sup>:

| Log-likelihood | Threshold fold | Exact rule | Difference |
| --- | ---: | ---: | ---: |
| `log_lh_seq` | -341,165.61 | -341,139.05 | +26.6 |
| `log_lh_time` | -48,870.21 | -48,757.55 | +112.7 |
| `log_lh_coal` | -24,152.82 | -24,116.64 | +36.2 |

The exact rule moves 1,147 node dates by more than 1e-9 years (99th percentile 0.0012 years, maximum 0.008 years):

| Shift of the node date | Nodes |
| --- | ---: |
| None | 5,249 |
| Less than 15 minutes | 428 |
| 15 minutes to 2.4 hours | 351 |
| 2.4 hours to 1 day | 352 |
| 1 to 3 days | 16 |

- **Direction**: 735 nodes move later (mean +0.18 days) and 412 earlier (mean -0.07 days). The 16 largest shifts are all later (1.7 to 3.0 days): the root, 2 tips without an exact date, and 13 internal nodes. In this build, the fold therefore places dates systematically earlier than the exact rule
- **Root**: the root date moves 2.54 days later. Its <a id="gloss-use-5"></a>90% interval <sup>[5](#gloss-5)</sup> with the threshold fold is 0.28 days wide and does not contain the exact date. One more node next to the root has the same result; the other 1,145 moved nodes keep the exact date inside their 90% interval
- **Scale**: the clock rate gives $5.7 \cdot 10^{-5} \cdot 197{,}209 = 11.2$ substitutions per genome per year, one substitution about every 33 days. The largest shift, 2.95 days, is about one tenth of that time
- **Path to the dates**: the mutation lengths of the output are identical in both runs. The fold changes the dates through the branch-length distributions that the time inference computes from the sequence partition

In one forward pass, the fold removes about 5.0M explicit posterior entries, summed over the nodes. The exact rule keeps them; about 4.4k of them differ from the shared profile by 1e-9 or more (2 by more than 0.1, about 2.4k by more than 1e-6), about 1.0M by 1e-12 to 1e-9, 0.95M by 1e-16 to 1e-12, and 3.06M by less than 1e-16.

## Cost of the bitwise fold at a root polytomy

At the root, `fn is_shared_profile()` folds a site where a child disagrees only when its profile equals the shared profile of its state bit for bit. The site lacks the term of a child without data there, or carries the term of another state, so its subleading probabilities differ from those of the shared profile until both underflow to a one-hot profile. Below the root, the forward pass folds the same sites at $1 - 10^{-4}$. The branch-length contribution of an edge groups its sites by the bit pattern of their coefficients (`struct GroupedSiteContributions` in [packages/treetime/src/partition/optimize/sparse.rs](../../packages/treetime/src/partition/optimize/sparse.rs)), so an explicit site costs an extra entry only when its subleading probabilities reach the last bits of the coefficients.

`test_root_polytomy_child_contributions_evaluate_the_sites_of_the_same_polytomy_below_the_root` in [packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs](../../packages/treetime/src/partition/marginal/sparse/__tests__/test_root_polytomy.rs) compares the contribution entries of the child edges of a star at the root and one level below it. Each tip except `c0` has 2 `N` runs:

| Case | Children, branch length, sites | Explicit root sites, star at the root | Explicit sites of `P`, star below the root | Child contribution entries, at the root | Child contribution entries, below the root |
| --- | --- | ---: | ---: | ---: | ---: |
| `short_branches` | 12, 0.01, 60 | 58 | 1 | 114 | 67 |
| `medium_star` | 50, 0.01, 60 | 60 | 0 | 240 | 240 |
| `large_star` | 150, 1e-4, 120 | 0 | 0 | 732 | 732 |

- **12 children**: the subleading probabilities of the root sites are 1e-25 to 1e-15 and change the coefficients, so the contributions at the root hold 1.7 times the entries below the root. The case is ignored
- **50 children**: the root keeps every site explicit, but the subleading probabilities are 1e-112 to 1e-70, the coefficients are bit-identical, and the grouping merges them
- **150 children**: the profiles underflow to one-hot profiles and fold at both placements, as on the mpox B.1 build with its 657 root children

The class-keyed shared profiles of [kb/proposals/sparse-shared-profiles-for-sites-without-data.md](../proposals/sparse-shared-profiles-for-sites-without-data.md) remove both effects: sites with the same state and the same pattern of disagreeing children share one exact profile, so the root needs no bitwise equality and the forward pass needs no threshold.

## Size of the differences

TreeTime log-likelihoods use the natural logarithm, so their unit is the nat. A difference of $d$ nats between two log-likelihoods means that one likelihood is $e^d$ times the other. Log-likelihoods add up over sites and branches, so small differences per site give a large total.

| Difference | Likelihood ratio | Example |
| ---: | ---: | --- |
| 1e-6 nats | 1.000001 | Tolerance of the sparse-against-dense tests |
| 1e-4 nats | 1.0001 | Size of the fold threshold: $\ln(1 - 10^{-4}) \approx -10^{-4}$ |
| 0.04 nats | 1.04 | One coin toss that lands heads: a coin with 52% heads explains it 4% better than a fair coin, $\ln(0.52 / 0.50) = 0.039$. The largest edge objective difference of `long_branches_below_root` (0.043 nats) has this size |
| 0.69 nats | 2 | One bit: a fair coin toss has likelihood 1/2, and $\ln 2 = 0.69$ |
| 1.92 nats | 6.8 | Drop from the maximum at the edge of a 95% confidence interval for one parameter: half of 3.84, the 95% quantile of $\chi^2_1$ <a id="cite-1"></a>[Wilks 1938](https://doi.org/10.1214/aoms/1177732360) [[1](#ref-1)] |
| 26.6 nats | 3.6e11 | Change of `log_lh_seq` of the mpox B.1 build above |

### Why a fold below 1e-4 gives errors above 1e-6

The fold replaces the exact profile of a site by the shared profile of its state. The absolute change of each probability is below $10^{-4}$, but the relative change of the small probabilities can be large. At node `P`, position 21, of `long_branches_below_root`, the probability of the other states is $8 \cdot 10^{-5}$ in dense and $5 \cdot 10^{-6}$ in sparse, a factor of 16. The outside message of every child of `P` carries this change into its edge objective.

The size of the change depends on the branch lengths. With long branches, each child carries little information about the state of its parent, so the exact posterior without the information of one child differs more from the shared profile, which assumes that every child agrees. With branch length 0.01, the same tree gives differences far below 1e-6.

### Assessment

- **Dates of the mpox B.1 build**: small. 99.7% of the nodes move by less than one day, and 1,145 of the 1,147 moved nodes keep the exact date inside their 90% interval
- **One branch**: small on the scale of statistical support. 0.043 nats is about 2% of the 1.92 nats that bound a 95% confidence interval of one parameter
- **Contract of the sparse representation**: a defect. The difference is about 40,000 times the 1e-6 tolerance, it grows with the branch length and with the number of children that disagree, and no analysis bounds it
- **Comparisons of total log-likelihoods**: misleading. Differences of 26.6 to 112.7 nats are large where model or run comparisons decide on a few nats
- **Root date**: the most affected output. It moves by 9 times the width of its 90% interval

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

The exact rule: the node's message to its parent enters the forward product under a fold rule that lets its explicit entries fold only when they equal the shared profile bit for bit, while the message from the parent keeps the threshold fold. The backward pass already applied the fold rule at the node with the references of all children, so the posterior of a node below the root then keeps the positions that the root keeps. Commit `7b4e0ba78` on branch `fix/sparse-posterior-below-root` implements it (not landed) and enables the ignored cases.

The exact rule doubles the run time of the production mpox B.1 build: 552 s against 270 s, one thread, both runs under the same concurrent load; on a quiet build host (load 0.5 at the start, 7.9 at the end), the approximate fold of commit `204b81185` takes 260.6 s, against 230.6 s for commit `2d813264e`. In one forward pass, the explicit posterior entries grow from 15.3k to 5.02M and the outside-message entries from 5.05M to 33.0M, although only 62k and 636k of them are bit-distinct. [kb/proposals/sparse-shared-profiles-for-sites-without-data.md](../proposals/sparse-shared-profiles-for-sites-without-data.md) describes a sharing mechanism that removes these copies without a tolerance.

> [!IMPORTANT]
> **Decision required.** `rust` keeps the threshold fold below the root until the decision on the proposal. Options: take the exact rule with single-site entries (sparse matches dense, the production build takes twice as long), take it together with the class-keyed shared profiles of the proposal (exact, not implemented), or keep the fold.

## Related issues

- [H-ancestral-sparse-dense-internal-residues-diverge.md](H-ancestral-sparse-dense-internal-residues-diverge.md): names the fold of resolved positions as one lead for the internal-node residues that differ between the representations
- [M-sparse-dense-likelihood-residual-without-missing-data.md](M-sparse-dense-likelihood-residual-without-missing-data.md): a likelihood difference without missing data, where substitution sites fold below the root

## Glossary

1. <a id="gloss-1"></a> **Shared profile.** The probability vector that sparse mode stores once per state for all sites of a node that it does not store explicitly. It assumes that every child of the node agrees with that state. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Outside message.** The probability of the parent states of a node given all data outside the subtree of the node. Sparse mode computes it as the parent posterior divided by the node's backward message, the subtree likelihood of the pruning algorithm <a id="cite-2"></a>[Felsenstein 1981](https://doi.org/10.1007/bf01734359) [[2](#ref-2)]. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Edge objective.** The log-likelihood of the data as a function of the length of one branch, with the messages at both ends of the branch held fixed. Branch-length optimization maximizes it; the test compares its change at three branch lengths and its optimum. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Nat.** The unit of a logarithm to base $e$. A log-likelihood difference of $d$ nats is a likelihood ratio of $e^d$; 1 nat is a ratio of 2.72, and $\ln 2 = 0.69$ nats is one bit. [↩](#gloss-use-4)
5. <a id="gloss-5"></a> **90% interval.** The date range that holds 90% of the posterior probability of a node date, the lower and upper bounds of `timetree.dates.tsv`. [↩](#gloss-use-5)

## References

1. <a id="ref-1"></a> Wilks, S. S. 1938. "The Large-Sample Distribution of the Likelihood Ratio for Testing Composite Hypotheses." _The Annals of Mathematical Statistics_ 9 (1): 60-62. https://doi.org/10.1214/aoms/1177732360 [↩](#cite-1)
2. <a id="ref-2"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17 (6): 368-376. https://doi.org/10.1007/bf01734359 [↩](#cite-2)
