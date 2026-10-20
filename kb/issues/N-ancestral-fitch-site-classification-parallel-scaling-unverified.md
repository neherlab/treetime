# Parallel scaling of the Fitch setup is unverified

> [!IMPORTANT]
> **Investigation required.** Source inspection does not establish that complete-command scaling recovered, or that the parallel leaf conversion costs time at one worker. Collect the measurements listed under "Evidence required".

The indexed Fitch implementation previously classified sites serially before parallel reconstruction. Representative benchmarks at commit `20164657c95d853f15e1a5dc7d8955c6134798ce` showed slower eight- and sixteen-worker ancestral, optimize, and timetree commands despite improving isolated single-worker work.

`FitchSiteIndex` and its serial classification pass no longer exist. `fn attach_seqs_to_graph()` now looks up the leaf sequences serially and converts each leaf with `FitchNodeData::new()` in a Rayon `into_par_iter()` [packages/treetime/src/partition/fitch/passes.rs#L53-L88](../../packages/treetime/src/partition/fitch/passes.rs#L53-L88). The Fitch backward and forward passes run on the tree scans: `fn fitch_backward()` through `GraphSchedule::scan_up()` and `fn fitch_forward()` through `GraphSchedule::scan_down()` [packages/treetime/src/partition/fitch/passes.rs#L90-L218](../../packages/treetime/src/partition/fitch/passes.rs#L90-L218).

A single earlier measurement on `mpox/clade-ii/2000` found that the parallel leaf setup increased the ancestral wall time at `-j 1` from 5.732 s to 6.057 s (5.7%), while it improved the wall time at 2 to 16 workers; optimize and timetree showed no comparable one-worker regression. That run used one warmup and three measured runs per binary and configuration, on code that also parallelized the leaf lookup, so it does not establish a systematic regression of the current code.

## Evidence required

- Run representative complete ancestral, optimize, and timetree commands at one, eight, and sixteen workers from one revision, with `just profile` or the committed benchmark drivers
- Repeat paired one-worker measurements of sparse `ancestral` with and without the parallel leaf conversion, on mpox and on at least one dataset with a different compression ratio, and check whether the difference exceeds the run-to-run variation
- Compare outputs for semantic equivalence and deterministic ordering
- Attribute any remaining scaling loss, and compare any proposed one-worker path for output and error-state equivalence, before implementation starts

## Related issues

- [N-ancestral-parallel-sparse-leaf-validation-coverage.md](N-ancestral-parallel-sparse-leaf-validation-coverage.md)
