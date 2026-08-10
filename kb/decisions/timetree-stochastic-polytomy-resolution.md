# Timetree uses corrected stochastic polytomy resolution

v1 resolves temporal polytomies with one mutation-conditioned stochastic event process. It does
not provide v0's greedy pairwise method or an exact v0 random stream.

**Type**: Intentional v0 divergence with three v0 erratum corrections.

**v1 location**:
[`packages/treetime/src/timetree/optimization/polytomy/`](../../packages/treetime/src/timetree/optimization/polytomy/),
with merger-rate schedules assembled in
[`packages/treetime/src/timetree/pipeline.rs`](../../packages/treetime/src/timetree/pipeline.rs).

## Event model

The sweep moves from the most recent child toward the polytomy parent. A child becomes live when
the sweep reaches its calendar time. A live lineage can merge only after all substitutions mapped
to its incoming branch have received mutation events.

Let $M$ be the total remaining substitution count and let $R$ be the live mutation-free lineages.
The event rates are

$$R_{\mathrm{mut}}=\mu M, \qquad R_{\mathrm{coal}}=\max(0, |R|-1)\kappa(t),$$

where $\mu$ is the whole-alignment mutation rate and $\kappa(t)$ is the per-branch merger rate from
the active coalescent model. A mutation event selects a lineage in exact proportion to its integer
remaining count. A coalescent event selects two lineages uniformly from $R$.

Reconstructed edge substitutions are the primary count source. When they are unavailable, the
count is `round(mutation_length * alignment_length)`.

## Time-varying hazard

The sampler draws one unit-exponential hazard threshold for each stochastic event. It integrates
that threshold across every interval where the live-lineage set and $\kappa(t)$ are constant.
Lineage arrivals, merger-rate breakpoints, and the parent time are deterministic boundaries.

Crossing a lineage arrival or rate breakpoint carries the remaining hazard into the next
interval. Crossing the parent ends the sweep without an event. The remaining lineages stay as
children of the polytomy parent, so a short time window can retain a residual polytomy.

The active coalescent model supplies $\kappa(t)$. A skyline fit supplies its exact step schedule.
When no coalescent prior is requested, the pipeline estimates a constant $T_c$ from the tree and
uses its constant merger-rate schedule. This rate ownership is also defined in
[timetree-frozen-lineage-counts-for-coalescent-prior.md](timetree-frozen-lineage-counts-for-coalescent-prior.md).

## Random generator and command behavior

v1 uses the project Rust random generator. `--seed` makes the sampled topology reproducible; when
the user omits it, the generated seed is logged. Exact seeded parity with NumPy is not a target
because the corrected event process has a different distribution and draw order.

The timetree command has one resolution strategy. It does not reproduce v0's
`--stochastic-resolve` and `--greedy-resolve` method selection.

This event process samples plausible mutation-conditioned histories. It is not a calibrated
posterior over binary and unresolved topologies. Such a posterior requires a separate approved
model and validation contract.

## Correctness contract

- Input times and rates are finite; rates are non-negative and do not overflow.
- No merger is at or before the parent time or more recent than either child.
- Every lineage is referenced exactly once in the returned forest.
- A merger references only original children or earlier mergers.
- Graph mutation starts only after the complete plan passes validation.
- Reparenting preserves each existing child edge key and payload.
- The same input and random-generator state produce the same merger plan.

Unit tests cover event boundaries, zero rates, invalid rates, parent bounds, exact mutation weights,
calendar translation, plan validation, and graph application. Property tests cover plan and graph
invariants over arbitrary `u64` seeds. The first-merger waiting-time test checks the analytical
mean $1/((k-1)\kappa)$ for constant $\kappa$.

## v0 corrections

- [kb/v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md](../v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md)
- [kb/v0-errata/timetree-stochastic-resolve-event-past-parent.md](../v0-errata/timetree-stochastic-resolve-event-past-parent.md)
- [kb/v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md](../v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md)

## Related

- [kb/proposals/timetree-stochastic-polytomy-resolution.md](../proposals/timetree-stochastic-polytomy-resolution.md)
- [kb/algo/timetree.md](../algo/timetree.md)
- [kb/features/timetree.md](../features/timetree.md)
