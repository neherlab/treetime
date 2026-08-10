# generate_subtree applies events drawn past the parent node, producing negative branch lengths

## v0 location

`TreeTime.generate_subtree()` (`#TreeTime`, `#generate_subtree`) [packages/legacy/treetime/treetime/treetime.py#L872-L1010](../../packages/legacy/treetime/treetime/treetime.py#L872-L1010)

## Erratum

The sweep is bounded by the parent node's time: it may only place coalescent events in the
window between the most recent child and the parent. The bound is enforced by the loop
condition at [packages/legacy/treetime/treetime/treetime.py#L906](../../packages/legacy/treetime/treetime/treetime.py#L906):

```python
while len(branches_alive) + len(branches_to_come) > 2 and t < tmax:
```

but that test runs at the *top* of the iteration. Within the body, the drawn waiting time is
applied and the event is committed with no further check
([L928-L971](../../packages/legacy/treetime/treetime/treetime.py#L928-L971)):

```python
total_rate_inv = 1.0 / (total_mut_rate + total_coalescent_rate)
dt = exp_dis(total_rate_inv)
t += dt
...
    new_node.time_before_present = t
```

Nothing constrains `t` to remain below `tmax`. The exponential draw is unbounded, so on the
final iteration of any sweep `t` overshoots the parent with probability
$e^{-r\,(t_{\max} - t)}$, and a coalescent node is created at a time *older than the parent
node* before the loop condition notices.

The consequence surfaces in the epilogue at
[packages/legacy/treetime/treetime/treetime.py#L999-L1002](../../packages/legacy/treetime/treetime/treetime.py#L999-L1002):

```python
for b in branches_alive + branches_to_come:
    b.branch_length = tmax - b.time_before_present
```

For the overshooting node, `b.time_before_present > tmax`, so `branch_length` is **negative**.
The node is attached to the parent as a child sitting further in the past than its own parent,
inverting the time ordering of the tree at that edge.

The same overshoot also affects mutation events, though harmlessly: a mutation deducted at
`t > tmax` is a mutation placed outside the available window, which biases the count of
lineages that become eligible to coalesce but leaves no invalid value behind.

## Correct formulation

The overshoot means no event occurred within the window. The candidate time must be tested
against the bound before the event is committed; when it crosses, the sweep terminates and
the surviving lineages attach to the parent as a residual polytomy -- which is the same
outcome the loop condition produces one iteration later, minus the invalid node.

## Evidence

- The loop condition at [L906](../../packages/legacy/treetime/treetime/treetime.py#L906) establishes `t < tmax` as an invariant of the algorithm, but the body at [L930](../../packages/legacy/treetime/treetime/treetime.py#L930) can violate it before the next test
- The epilogue at [L1001](../../packages/legacy/treetime/treetime/treetime.py#L1001) computes `tmax - b.time_before_present` without a floor, so it assumes the invariant holds
- The adjacent greedy path handles the same constraint correctly: `_poly`'s merge-time optimization is bracketed to `[max_pos, min_pos]` derived from the parent and child times ([packages/legacy/treetime/treetime/treetime.py#L744-L762](../../packages/legacy/treetime/treetime/treetime.py#L744-L762)), so no merged node can be placed outside the window
- `branch_length` is used downstream as a non-negative quantity throughout `TreeTime`; no caller handles a negative value

## v0 impact

- One coalescent node per affected sweep receives a negative branch length and a time older
  than its parent
- The affected node's subtree is dated inconsistently until the next timetree inference pass
  reassigns node times
- Probability of occurrence per sweep is $e^{-r\,(t_{\max}-t)}$ evaluated at the last event,
  so it rises as the remaining window shrinks -- i.e. it is most likely exactly when the sweep
  is about to terminate normally

## v1 status

Implemented with the correction in [`simulate_subtree()`](../../packages/treetime/src/timetree/optimization/polytomy/sweep.rs). The parent time is an explicit deterministic hazard boundary. Crossing it ends the sweep without an event, and [`apply_plan()`](../../packages/treetime/src/timetree/optimization/polytomy/apply.rs) rejects any merger at or before the parent before graph mutation.
