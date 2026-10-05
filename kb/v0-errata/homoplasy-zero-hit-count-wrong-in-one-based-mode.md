# Homoplasy zero-hit site count is wrong when positions count from 1

## v0 location

`scan_homoplasies()` collects the branches hit at each position in `positions`, a `defaultdict(list)` keyed by `pos + offset`, where `offset` is 1 unless `--zero-based` is given [`packages/legacy/treetime/treetime/wrappers.py#L133-L156`](../../packages/legacy/treetime/treetime/wrappers.py#L133-L156). The per-taxon loop then reads `len(positions[pos])` with the zero-based key [`packages/legacy/treetime/treetime/wrappers.py#L158-L163`](../../packages/legacy/treetime/treetime/wrappers.py#L158-L163).

## Erratum

In the default one-based mode, the read `positions[pos]` of each terminal mutation at a site hit two or more times inserts an empty list under the zero-based key `pos` when that key was not hit, because `positions` is a `defaultdict`. `np.bincount` then counts these empty entries as sites with zero hits, and the next line overwrites the zero-hit bin with `L - np.sum(multiplicities_positions)` [`packages/legacy/treetime/treetime/wrappers.py#L183-L184`](../../packages/legacy/treetime/treetime/wrappers.py#L183-L184), which now subtracts the inserted keys as well. The number of sites with zero hits drops by the number of inserted keys, and the log-likelihood difference to the Poisson distribution changes with it.

## Evidence

- The two indexing modes describe the same reconstruction and differ only in the printed positions, yet they print different site counts. On `data/zika/86` with `--gtr jc`, the default run prints `10175 were hit 0 times` and a log-likelihood difference of `-6.666e+01`; the `--zero-based` run prints `10241 were hit 0 times` and `-7.077e+01`
- The rows must sum to the number of sites: $10807 - 480 - 69 - 12 - 5 = 10241$. The one-based rows sum to 10741
- The adjacent condition of the same loop uses the correct key: `pos + offset in positions and len(positions[pos + offset]) > 1` [`packages/legacy/treetime/treetime/wrappers.py#L161`](../../packages/legacy/treetime/treetime/wrappers.py#L161)
- The v0 tutorial [`packages/legacy/treetime/docs/source/tutorials/homoplasy.rst#L37-L41`](../../packages/legacy/treetime/docs/source/tutorials/homoplasy.rst#L37-L41) prints the wrong value 10175

## v0 impact

Every default `treetime homoplasy` run reports too few sites without mutations and a wrong log-likelihood difference. The size of the error depends on how many taxa carry mutations at hit positions.

## v1 status

v1 computes the zero-hit count as the number of sites minus the sites with at least one substitution, independent of the indexing mode (`fn site_histogram()` in [`packages/treetime/src/homoplasy/site_hits.rs`](../../packages/treetime/src/homoplasy/site_hits.rs)). On `data/zika/86`, v1 prints 10241 in both modes, as recorded in [kb/decisions/homoplasy-mutation-mapping-and-counting.md](../decisions/homoplasy-mutation-mapping-and-counting.md).
