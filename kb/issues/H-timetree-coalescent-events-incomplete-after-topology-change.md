# Coalescent model may be built before node times are recomputed after a topology change

`run_timetree` builds the coalescent model before its own backward/forward passes recompute node times: `compute_coalescent_model` runs at [packages/treetime/src/timetree/inference/runner.rs#L53-L55](../../packages/treetime/src/timetree/inference/runner.rs#L53-L55), before `propagate_distributions_backward`/`propagate_distributions_forward` at [runner.rs#L60-L64](../../packages/treetime/src/timetree/inference/runner.rs#L60-L64). The coalescent model therefore relies on node times left over from a previous pass. Any topology mutation that clears internal `time_distribution`s and is not followed by a coalescent-free time pass leaves the model built from an incomplete tree.

## Fixed

- The post-ancestral reroot instance is fixed: the pipeline now runs a coalescent-free `run_timetree` right after the reroot, before the optimization loop, so round 1's coalescent build sees a complete tree ([packages/treetime/src/timetree/pipeline.rs#L281-L293](../../packages/treetime/src/timetree/pipeline.rs#L281-L293)).
- Event collection no longer silently drops nodes without an inferred time. `collect_tree_events` returns a contextual error naming the offending node ([packages/treetime/src/coalescent/events.rs#L27-L54](../../packages/treetime/src/coalescent/events.rs#L27-L54)), so any remaining instance fails loudly instead of producing a wrong coalescent prior.

## Remaining

The resolve-polytomies refinement branch clears internal `time_distribution`s via `prepare_tree_after_topology_change` ([packages/treetime/src/timetree/optimization/polytomy.rs#L560](../../packages/treetime/src/timetree/optimization/polytomy.rs#L560)) and then calls `run_timetree` directly ([packages/treetime/src/timetree/refinement.rs#L74-L81](../../packages/treetime/src/timetree/refinement.rs#L74-L81)). That `run_timetree` builds the coalescent model before recomputing node times, so a dataset that actually resolves polytomies while a coalescent prior is active can still fail. This has not been reproduced (the tested datasets did not resolve polytomies under the coalescent modes), and it now surfaces as the explicit error from `collect_tree_events` rather than silent incompleteness.

## Potential solution

Recompute node times immediately after any topology mutation that clears them and before the coalescent model is consumed, or restructure `run_timetree` so the coalescent model is built from node times established for the current topology rather than a prior pass.

## Validation

- Coalescent modes on inversion-prone datasets (zika/20, flu/h3n2/200) no longer abort with "Lineage count must end at zero".
- Inject one missing internal time and require a contextual error naming the node ([packages/treetime/src/coalescent/**tests**/test_events.rs](../../packages/treetime/src/coalescent/__tests__/test_events.rs)).
- Exercise a resolve-polytomies run that actually introduces nodes while a coalescent prior is active, and assert the coalescent model is built from complete node times.
