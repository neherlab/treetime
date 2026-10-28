# Shared profiles for sparse sites without data below a child

Sparse marginal inference stores a message as one profile per variable site (`variable`) and one profile per state with a count for every other site (`fixed`). A site where a child subtree has only gaps or `N` cannot use the `fixed` profile of its state, because the child contributes no factor there. The sparse passes therefore keep each such site as a single-site entry in every message above that subtree. Alignments with long `N` runs produce thousands of such entries per node, and sparse `optimize` becomes slower and needs more memory.

This proposal stores these sites in shared profiles again, one per pattern of children without data, and records the measurements of both designs.

## Current behavior

- The backward pass adds every site in a child's no-data range that lies outside the node's own gaps and unknown sites as a candidate site of the node (`fn combine_child_messages()` [packages/treetime/src/partition/marginal/sparse/backward.rs#L133-L206](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L133-L206)). Each candidate site gets its own profile in `variable`
- The forward pass walks the same sites one by one to remove them from the per-state counts of the outside message (`fn compute_msg_to_child()` [packages/treetime/src/partition/marginal/sparse/forward.rs#L275](../../packages/treetime/src/partition/marginal/sparse/forward.rs#L275))
- The no-data range of each child travels with the backward edge message (`child_no_data` of `struct SparseEdgeBackward` [packages/treetime/src/partition/storage/sparse.rs#L130](../../packages/treetime/src/partition/storage/sparse.rs#L130))

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

## Trade-offs

- **Shared profiles**: time and memory stay close to the code before the fix. In exchange, every reader of sparse messages handles a third storage kind, and the grouping logic adds about 600 lines
- **Single-site entries**: one storage kind and the smallest code. In exchange, time and memory grow with the number of sites without data below a child; on the largest example datasets, runs take seconds longer and need up to four times the memory
- **Single-site entries with a cheaper representation**: for example, store the no-data ranges once per node and evaluate their contribution in closed form per state and range length. This keeps one storage kind in the messages but needs its own derivation for the forward pass and the objective, and nobody has implemented or measured it

> [!IMPORTANT]
> **Decision required.** Keep single-site entries, reintroduce the shared profiles of commit `107ec665b`, or design a cheaper single-site representation. The measurements above are the evidence; no option is approved.
