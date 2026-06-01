# --divergence-units=mutations not implemented for node data divergence

Augur's refine accepts `--divergence-units {mutations-per-site, mutations}`, controlling whether the divergence that `augur export v2` accumulates into `node_attrs.div` is reported in substitutions per site or in absolute mutation counts (per-site values scaled by alignment length).

Timetree node data JSON emits `mutation_length`, and optimize node data JSON emits `branch_length`, only in substitutions per site (the `mutations-per-site` convention). There is no flag to request absolute mutation counts, so a consumer expecting `mutations` units receives per-site values without conversion. In augur's non-timetree refine path the `mutations` mode overwrites `branch_length` with an integer count (`refine.py`); v1's optimize writer does not implement this.

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
