# generate_subtree skips the interval between a branch arrival and the overshooting draw

## v0 location

`TreeTime.generate_subtree()` (`#TreeTime`, `#generate_subtree`) [packages/legacy/treetime/treetime/treetime.py#L872-L1010](../../packages/legacy/treetime/treetime/treetime.py#L872-L1010)

## Erratum

The sweep advances backwards in time through a set of branches that become available
progressively: `branches_to_come` are ordered by node time, and each becomes eligible to
participate once the sweep reaches it. Because the event rate depends on how many branches
are alive, the rate is piecewise constant with breakpoints at the arrival times, and a
waiting time drawn under one rate is only valid up to the next breakpoint.

v0 handles the breakpoint at
[packages/legacy/treetime/treetime/treetime.py#L929-L936](../../packages/legacy/treetime/treetime/treetime.py#L929-L936):

```python
dt = exp_dis(total_rate_inv)
t += dt
# if the time advanced past the next branch in the branches_to_come list
# add this branch to branches alive and re-renter the loop
if len(branches_to_come) and t > branches_to_come[0].time_before_present:
    while len(branches_to_come) and t > branches_to_come[0].time_before_present:
        branches_alive.append(branches_to_come.pop(0))
```

Detecting the crossing and discarding the event is correct; the exponential distribution is
memoryless, so redrawing under the new rate is valid. But the resumption point is wrong.
`t` has already been advanced to `t + dt`, which lies *beyond* the arrival time
$t_a$, and the loop re-enters from there.

The interval $[t_a,\ t + dt]$ is therefore never evaluated at any rate. It is skipped
entirely: no event can be drawn in it, even though the higher post-arrival rate applies
throughout. The correct resumption point is $t_a$ itself.

Because arrivals only ever *increase* the number of alive branches, the skipped interval is
always one where the rate is higher than the rate used to draw across it. The bias is
therefore systematic, not symmetric: waiting times are overestimated and events are pushed
further into the past than the model specifies.

The inner `while` compounds it. If the single draw overshoots several arrivals, all of them
are consumed at once and the whole span from the first arrival to `t + dt` is skipped, with
the rate rising at each arrival passed.

## Correct formulation

On detecting a crossing, set `t` to the arrival time of the branch being admitted, admit that
branch (and any tied at the same time), and redraw. Each interval between consecutive
breakpoints is then evaluated exactly once, at its own rate.

## Evidence

- The comment at [L931-L932](../../packages/legacy/treetime/treetime/treetime.py#L931-L932) states the intent as "add this branch to branches alive and re-renter the loop", i.e. treat the crossing as a rate change rather than as elapsed time -- which is what resuming at $t_a$ accomplishes and resuming at `t + dt` does not
- The event is correctly discarded rather than applied, showing the author recognised the drawn time as invalid beyond the breakpoint; the position update is not rolled back to match
- `dummy_coalescent_rate` at [L896](../../packages/legacy/treetime/treetime/treetime.py#L896) is calibrated as $2/(t_{\max} - t)$ so that the lineages "typically coalesce" within the available window ([L894-L895](../../packages/legacy/treetime/treetime/treetime.py#L894-L895)); systematically overestimating waiting times works against that stated calibration
- The rate recomputation at the top of the loop ([L908-L920](../../packages/legacy/treetime/treetime/treetime.py#L908-L920)) is unconditional, confirming the design treats the rate as piecewise constant between arrivals

## v0 impact

- Coalescence and mutation events are placed systematically further into the past than the
  specified process, by an amount that grows with the number of children and the spread of
  their node times
- Fewer events fit inside the window before the parent bound is reached, so polytomies are
  left less resolved than the model implies
- Largest polytomies are worst affected: more children means more arrivals, hence more skipped
  intervals per sweep

## v1 status

Implemented with the correction in [`simulate_subtree()`](../../packages/treetime/src/timetree/optimization/polytomy/sweep.rs). One unit-exponential hazard threshold is carried across every lineage-arrival and merger-rate boundary. The sweep evaluates each interval under its own constant rate and resumes exactly at the crossed boundary.
