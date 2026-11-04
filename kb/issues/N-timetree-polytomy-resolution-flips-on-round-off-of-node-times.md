# Polytomy resolution flips on round-off changes of node times

## Problem

The stochastic polytomy resolution of `timetree` ([kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md)) samples the history of each polytomy from the seeded random stream of the run, with the node times of the first iteration as input. The sample is a discrete outcome of continuous inputs: the node times set the time window and the order in which children become live, and the comparisons of sampled event times against them decide which pairs merge. A round-off change of a node time can therefore change which pairs merge, and with them the topology of the output and the dates of the nodes nearby. Two builds of the same commit with different floating-point round-off (another compiler, CPU, or FFT plan) do not give the same tree for the same `--seed`.

`fn resolve_polytomies()` ([packages/treetime/src/timetree/optimization/polytomy/resolve.rs](../../packages/treetime/src/timetree/optimization/polytomy/resolve.rs)) passes one random stream through all polytomies in turn. When the flipped polytomy draws a different number of random numbers, every polytomy after it in that order reads a shifted stream and can resolve differently, too.

## Evidence

On the Nextstrain lassa L build, two runs whose node times differ by round-off (second skyline log-likelihood -6930.9874 in one run and -6930.9868 in the other) resolve every polytomy alike except one of the 41. The four-child polytomy at `NODE_0000423` merges `MK107847` and `MK107849` in one run and merges no pair in the other, so the 2-tip clade `NODE_0000752` exists only in the first run, and the parent of the polytomy moves by 0.080 years. Source: "Exact tails, production runs" in [kb/decisions/distribution-convolution-tails-by-exponential-tilting.md](../decisions/distribution-convolution-tails-by-exponential-tilting.md).

## v0

v0 has the same sensitivity in its stochastic method (`--stochastic-resolve`, `TreeTime.generate_subtree()` [packages/legacy/treetime/treetime/treetime.py#L872-L1010](../../packages/legacy/treetime/treetime/treetime.py#L872-L1010)): it draws exponential waiting times and the merging pairs from the one generator `self.rng`, starts from the `time_before_present` of the children, and compares the sampled time with the next child time, so a round-off change of a child time can change the outcome and shift the stream for the later polytomies. The greedy method of v0 (`TreeTime._poly()`, the default without `--stochastic-resolve`) draws no random numbers, but it merges the pair with the largest cost gain while that gain exceeds `resolution_threshold = 0.05`, so a near tie or a gain near the threshold flips in the same way.

## Fix direction

Give each polytomy its own random stream, derived from the run seed and a stable key of the polytomy (for example the name of the polytomy node), so that a flip in one polytomy no longer shifts the samples of the others. The flip of a single polytomy remains: it is the expected behavior of a sample whose inputs changed, of the same kind as a change of `--seed`.

## Validation

A test with two polytomies whose node times differ by one ulp at the first polytomy only shows that the second polytomy resolves alike in both runs.
