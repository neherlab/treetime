# Fitch site-classification parallel scaling is unverified

The indexed Fitch implementation previously classified sites serially before parallel reconstruction. Representative benchmarks at commit `20164657c95d853f15e1a5dc7d8955c6134798ce` showed slower eight- and sixteen-worker ancestral, optimize, and timetree commands despite improving isolated single-worker work.

Current [`attach_seqs_to_graph()`](../../packages/treetime/src/ancestral/fitch.rs) uses Rayon for leaf lookup and per-partition sparse node construction. Source inspection establishes that the original `FitchSiteIndex` serial classification pass is gone, but does not establish that site-level work is parallel or that current complete-command scaling recovered.

## Evidence required

- Re-run representative complete ancestral, optimize, and timetree commands at one, eight, and sixteen workers from one revision.
- Compare outputs for semantic equivalence and preserve deterministic ordering.
- Attribute any remaining scaling loss before creating an implementation ticket.

## Related issues

- [N-ancestral-parallel-sparse-leaf-single-thread-regression.md](N-ancestral-parallel-sparse-leaf-single-thread-regression.md)
- [M-benchmark-reports-mix-revisions-and-are-not-reproducible.md](M-benchmark-reports-mix-revisions-and-are-not-reproducible.md)
