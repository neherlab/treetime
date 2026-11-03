# Shared profiles for sparse sites without data below a child

Sparse marginal inference stores a message as one profile per variable site (`variable`) and one profile per state with a count for every other site (`fixed`). A site where a child subtree has only gaps or `N` cannot use the `fixed` profile of its state, because the child contributes no factor there. The sparse passes therefore keep each such site as a single-site entry in every message above that subtree. Alignments with long `N` runs produce thousands of such entries per node, and sparse `optimize` becomes slower and needs more memory.

This proposal stores these sites in <a id="gloss-use-1"></a>shared profiles <sup>[1](#gloss-1)</sup> again, one per pattern of children without data, and records the measurements of both designs. It also records the cost of exact posteriors below the root ([kb/issues/M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md](../issues/M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md)), which the same kind of sharing can remove.

## Current behavior

- The backward pass adds every site in a child's no-data range that lies outside the node's own gaps and unknown sites as a candidate site of the node (`fn combine_child_messages()` [packages/treetime/src/partition/marginal/sparse/backward.rs#L146-L219](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L146-L219)). Each candidate site gets its own profile in `variable`
- The forward pass walks the same sites one by one to remove them from the per-state counts of the <a id="gloss-use-2"></a>outside message <sup>[2](#gloss-2)</sup> (`fn compute_msg_to_child()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L285](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L285))
- The no-data range of each child travels with the backward edge message (`child_no_data` of `struct SparseEdgeBackward` [packages/treetime/src/partition/storage/sparse.rs#L142](../../packages/treetime/src/partition/storage/sparse.rs#L142))

The cost grows with the number of sites without data below a child, summed over all nodes, so it grows with tree size and with the length of the `N` and gap runs.

## Proposed design

Sites below a node that lack data in the same set of children, and that have the same Fitch state at the node, receive the same messages. Their profiles are therefore equal, and one profile with a count represents all of them, as `fixed` does for sites without variation.

- **Storage**: a third part of `SparseSeqDistribution` next to `variable` and `fixed`: a list of site ranges, each with a group number, and per group one profile and one count per state
- **Lookup**: the profile of a site and state comes from `variable` if the site is there, else from the group that covers the site, else from `fixed`
- **Building the groups at a node**: for each incoming message, mark each site range as covered by a group of that message, as without data below, or as plain. The tuple of these marks over all incoming messages is the group key. Ranges whose marks are all plain stay in `fixed`. For each remaining key, count the node's Fitch states over its sites, leaving out the candidate sites
- **Combining**: each group contributes one product per state, weighted by its count, to the log-likelihood and to the per-state counts, in the same way as the `fixed` profiles
- **Readers**: propagation along an edge, the outside message of the forward pass, the branch-length objective, and the readers of reconstructed sequences, mutations, and ambiguity codes each read the groups next to `variable` and `fixed`
- **Per-site rates**: with a rate per site, two sites never share a profile, so grouping stays off and these sites remain single-site entries

A full implementation exists in commit `107ec665b` (`perf(sparse): group sites without data below a child into shared profiles`), and commit `2a41a522e` removes it. It keeps the group storage and the group building in their own modules (`site_groups.rs`, `site_grouping.rs`, about 260 lines) and adds grouped branches to each reader listed above, about 600 lines in total.

## Measurements

### Method

- **Command**: `treetime optimize --representation=sparse --tree=<dataset>/tree.nwk --alignment=<dataset>/aln.fasta.xz --output-all=<dir>`, with default settings otherwise
- **Builds**: release profile (`just build-release treetime`) of three commits:
  - Before the fix: `a973a2f45` on `rust`. Its results are wrong at sites without data below a child
  - With grouping: `39a938bcb`
  - Without grouping: `2a41a522e`, the current state
- **Machine**: a Linux build host with 20 hardware threads, at a load average of 1.5 to 5 from other jobs
- **Quantities**: wall time, CPU time (user plus system), and peak resident memory from `/usr/bin/time -f "%e %U %S %M"`
- **Repetitions**: three alternating runs per binary for `sc2/2844` and `mpox/clade-ii/500` (the table gives the median), and one run per binary for the larger datasets

### Time and memory

| Dataset | Before the fix | With grouping | Without grouping |
|---|---|---|---|
| `dengue/2000` | 0.75 s, 5.33 s CPU, 200 MB | 1.19 s, 8.53 s CPU, 215 MB | 1.44 s, 11.89 s CPU, 281 MB |
| `sc2/2844` | 1.14 s, 3.45 s CPU, 220 MB | 1.66 s, 6.81 s CPU, 260 MB | 4.01 s, 27.93 s CPU, 568 MB |
| `sc2/5100` | 3.03 s, 8.30 s CPU, 469 MB | 5.79 s, 20.38 s CPU, 490 MB | 12.46 s, 118.30 s CPU, 1590 MB |
| `mpox/clade-ii/500` | 0.59 s, 2.47 s CPU, 227 MB | 0.93 s, 5.73 s CPU, 240 MB | 5.46 s, 60.48 s CPU, 964 MB |
| `mpox/clade-ii/2000` | 2.96 s, 12.51 s CPU, 863 MB | 4.20 s, 25.27 s CPU, 915 MB | 21.68 s, 287.79 s CPU, 3277 MB |

Wall time depends on the number of threads; the CPU time shows the total work. On a machine with 4 threads, `mpox/clade-ii/2000` without grouping takes over a minute.

### Earlier measurement during development

An earlier benchmark ran on the same host at a load average of 13 to 60, with default representation, five alternating runs per binary, and an ncov build input from the Nextstrain workflow (not in `data/`, run with `--model infer --no-indels`). These numbers are noisier and come from intermediate states of the code:

| Measurement | Before the fix | Without grouping | With grouping |
|---|---|---|---|
| `sc2/2844` wall time | 0.65 s to 0.74 s | 2.39 s | 0.87 s |
| `sc2/2844` CPU time | 2.12 s | | 3.45 s |
| `sc2/2844` peak memory | 216 MB to 225 MB | 595 MB | 255 MB to 264 MB |
| `sc2/2844` time per `optimize` iteration | 20 ms | 300 ms | |
| ncov wall time | 2.00 s to 2.06 s | 3.71 s | 2.16 s |
| ncov CPU time | 5.50 s | | 7.21 s |

### Results with and without grouping

Grouping sums the log-likelihood of n equal sites as one product times n instead of n products, so the rounding differs. `optimize` outputs of the two designs:

| Run | Result |
|---|---|
| `optimize` on `flu/h3n2/500` | identical files |
| `ancestral --method-anc=marginal --representation=sparse --impute-missing-data` on `mpox/clade-ii/500` | identical files |
| `prune --collapse-short=1e-6 --collapse-empty` on `mpox/clade-ii/500` | identical files |
| `optimize` on `sc2/2844` | same tree; 1371 branch lengths differ, by at most 2.8e-9 relative |
| `optimize` on `mpox/clade-ii/500` | same tree; 417 branch lengths differ, by at most 1.6e-9 relative |
| `optimize --no-indels` on `mpox/clade-ii/500` | same tree; 299 branch lengths differ, by at most 2.6e-6 relative |

The branch-length optimizer stops at its convergence tolerance, so a difference in the last digits of the objective moves the result within that tolerance. The sparse-against-dense tests of missing data pass at `1e-10` with both designs.

## Exact posteriors below the root

Below the root, the forward pass folds a site into the shared profile of its state when its posterior exceeds $1 - 10^{-4}$, also where the backward pass kept the site explicit because a child disagrees with the node: the child has no data there, or a substitution on the child edge gives it another state ([kb/issues/M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md](../issues/M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md)). The exact rule keeps such a site explicit in the posterior unless it equals the shared profile bit for bit: the node's message to its parent enters the forward product under a fold rule that lets its explicit entries fold only on exact equality (commit `7b4e0ba78` on branch `worktree/fix/sparse-posterior-below-root`, not landed). With it, sparse matches dense within 2e-12 on the root-polytomy test cases.

The exact rule keeps every one of these sites as a single-site entry in the posterior and in the outside message of every child. On a production mpox B.1 build (4,946 tips, 197,209 sites, a root polytomy of 657 children) it costs, in one forward pass:

| Quantity | Approximate fold | Exact rule |
|---|---:|---:|
| `timetree` wall time, one thread | 270 s | 552 s |
| Explicit posterior entries, summed over nodes | 15.3k | 5.02M |
| Explicit outside-message entries, summed over edges | 5.05M | 33.0M |
| Bit-distinct (state, profile) posterior entries, summed over nodes | | 62k |
| Bit-distinct (state, profile) outside entries, summed over edges | | 636k |

Almost all of the kept entries are bit-identical copies: sites with the same state and the same pattern of disagreeing children receive the same inputs in the same order. The difference to the shared profile is real but mostly tiny: about 4.4k entries differ from it by 1e-9 or more, and 3.06M by less than 1e-16. A sharing mechanism that represents equal inputs once removes the copies without a tolerance.

## Related techniques

- **Site repeats** (Kobert, Stamatakis, Flouri 2017, [doi:10.1093/sysbio/syw075](https://doi.org/10.1093/sysbio/syw075)): every directed partial likelihood vector carries a class number per site. At a node, the class of a site is a lookup of the pair of child classes, and each class is computed once. `corax_update_repeats()` in `.repos/coraxlib/src/corax/core/repeats.c` computes the lookup index `site_id_left[s] + site_id_right[s] * ids_left`, and `corax_default_enable_repeats()` turns repeats off for a node whose children have more than half as many classes as sites. RAxML-NG uses repeats by default (`opts.use_repeats = true` in `.repos/raxml-ng/src/CommandLineParser.cpp`) and turns them off for probabilistic alignments and for the tip-inner optimization. The method is exact, because equal inputs pass through the same operations. Its implementation in coraxlib handles two children per operation
- **CMAPLE genome lists** (`.repos/cmaple`): regions of reference, nucleotide, `N`, or explicit vectors. An `N` region of one child passes the other child's region through (`merge_N_RACGT_TwoLowers()` in `alignment/seqregions.cpp`), which is exact. Three steps are approximate: two identical reference regions merge into one with a log-likelihood term linear in the branch lengths (`merge_identicalRACGT_TwoLowers()`), an explicit vector becomes a single state when only one state exceeds `threshold_prob` (`SeqRegions::simplifyO()`, default `1e-8` in `utils/tools.cpp`), and region comparisons use tolerances. CMAPLE is in the same family as the fold at $1 - 10^{-4}$ and does not meet the exact target
- **Outside messages by division**: v1 divides the parent posterior by the child's message, and v0 does the same in log space (`node.marginal_outgroup_LH` in `packages/legacy/treetime/treetime/treeanc.py`). The messages from children are strictly positive for a branch length above the minimum of `fn fix_branch_length()`, and `fn divide_out_value()` guards underflow, so the division adds no cancellation. Prefix and suffix products avoid the division, but they keep the same set of explicit sites, because the outside message of a child depends on every sibling, and they compute the products of different children in different orders, so siblings no longer share bit-equal profiles. The division keeps the outside message of every agreeing child a function of the parent posterior and the child's message, which is the property that sharing needs

## Class-keyed shared profiles

A generalization of the groups of `107ec665b` covers the sites with substitutions too. A **class** is the state of a site together with the list of incoming messages whose value at the site differs from their shared profile, each with its own class. A message stores a table of classes with one profile per state and a list of position ranges that point to classes:

- The `fixed` profiles are the class of the empty list, and an explicit single site is a class with one position
- The groups of `107ec665b` are the classes whose list holds only children without data (labels `Fixed`, `Missing`, `Group(g)` of `fn segment_labels()` in `site_grouping.rs` of that commit)
- `combine_messages()` and `compute_msg_to_child()` step over the segments of the merged range boundaries of their inputs and compute each new class once, with counts per state for the log-likelihood and the branch-length weights. The branch-length objective already groups equal coefficients with a multiplicity (`struct GroupedSiteContributions` in [packages/treetime/src/partition/optimize/sparse.rs](../../packages/treetime/src/partition/optimize/sparse.rs))
- A site folds exactly when its class equals an existing class, so the fold needs no tolerance; per site, the profiles are those of the exact rule, and only the summation order of the log-likelihood changes
- With per-site rates, the rate is part of the class, so every site is its own class: no gain, and results stay exact

On the mpox B.1 build, the class count bounds the profile computations at about 62k per posterior pass and 636k for the outside messages, against 5.02M and 33.0M single-site entries.

> [!IMPORTANT]
> **Investigation required.** The cost of class-keyed profiles also grows with the number of range boundaries per node: the union of the `N`-run boundaries of all children, 657 at the root of mpox B.1. Nobody has counted these boundaries; count them on the mpox B.1 build before an implementation.

## Trade-offs

- **Shared profiles**: time and memory stay close to the code before the fix. In exchange, every reader of sparse messages handles a third storage kind, and the grouping logic adds about 600 lines
- **Single-site entries**: one storage kind and the smallest code. In exchange, time and memory grow with the number of sites without data below a child; on the largest example datasets, runs take seconds longer and need up to four times the memory
- **Exact posteriors below the root with single-site entries** (`7b4e0ba78`): sparse matches dense, and a production mpox B.1 build takes twice as long
- **Class-keyed shared profiles**: exact, and expected to remove both the cost of sites without data below a child and the cost of exact posteriors below the root. In exchange, every reader of sparse messages handles classes over ranges, which brings back a grouping mechanism of the size of `107ec665b` or larger; the range-boundary count is not measured
- **Single-site entries with a cheaper representation**: for example, store the no-data ranges once per node and evaluate their contribution in closed form per state and range length. This keeps one storage kind in the messages but needs its own derivation for the forward pass and the objective, and nobody has implemented or measured it

**Accuracy of the approximate fold.** The threshold fold that stays on `rust` until the decision changes sparse results below the root by up to 0.043 <a id="gloss-use-3"></a>nats <sup>[3](#gloss-3)</sup> per edge objective in the test polytomy with long branches. On a production mpox B.1 build, the exact rule raises `log_lh_seq` by 26.6 nats, moves 1,147 of 6,396 node dates (16 of them by 1 to 3 days), and moves the root date 2.54 days, out of the 90% interval of the threshold fold. Section "Size of the differences" of [kb/issues/M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md](../issues/M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md) explains the scale of these numbers.

> [!IMPORTANT]
> **Decision required.** The options are:
>
> - Keep single-site entries and the approximate fold below the root, the state of `rust`
> - Keep single-site entries and take the exact rule of `7b4e0ba78`
> - Implement class-keyed shared profiles with the exact rule, after the range-boundary count above
> - Reintroduce the shared profiles of commit `107ec665b` with the approximate fold below the root
> - Design a cheaper single-site representation
>
> The measurements above are the evidence; no option is approved. Until the decision, `rust` keeps the approximate fold below the root.

## Glossary

1. <a id="gloss-1"></a> **Shared profile.** A probability vector that sparse mode stores once, with a count, for many sites of a node instead of once per site. The `fixed` profiles hold one such vector per state for all sites that the node does not store explicitly, and assume that every child of the node agrees with that state. [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Outside message.** The probability of the parent states of a node given all data outside the subtree of the node. Sparse mode computes it as the parent posterior divided by the node's backward message. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Nat.** The unit of a logarithm to base $e$. A log-likelihood difference of $d$ nats is a likelihood ratio of $e^d$; 1 nat is a ratio of 2.72, and $\ln 2 = 0.69$ nats is one bit. [↩](#gloss-use-3)
