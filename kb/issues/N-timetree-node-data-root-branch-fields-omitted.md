# Root node omits placeholder branch-length fields in node data JSON

Augur's refine output writes branch-length fields (`branch_length`, `clock_length`, `mutation_length`) for every node, including the root, where they are placeholder values because the root has no parent edge.

Timetree node data JSON omits these fields on the root node. The root carries dates and confidence but no branch metrics. Downstream, `augur export v2` reads branch metrics from child edges when accumulating divergence and time, so the absent root values have no observed effect. The omission is a structural difference from augur's per-node field set.

## Related issues

- [Missing output files compared to v0](N-timetree-missing-output-files.md)
