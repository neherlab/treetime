# Marginal clock-rate probes report grid point limit errors without the explanation of time inference

When `timetree` estimates the clock rate with a coalescent prior, each refinement round steps the rate toward the maximum of the marginal likelihood of the node times. `fn update_clock_model()` calls `fn marginal_clock_rate()` when no clock rate is fixed and a coalescent prior exists ([packages/treetime/src/timetree/round.rs#L258-L270](../../packages/treetime/src/timetree/round.rs#L258-L270)), and `fn marginal_clock_rate()` calls `pub(crate) fn step_marginal_clock_rate()` ([packages/treetime/src/timetree/round.rs#L320-L351](../../packages/treetime/src/timetree/round.rs#L320-L351)). The step runs the backward pass three times, at the reference rate and at the rate times $e^{\pm h}$, through `pub(crate) fn marginal_log_lh_at_rate()` ([packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L17-L65](../../packages/treetime/src/timetree/optimization/marginal_clock_rate.rs#L17-L65)). Each probe rescales the branch-time distributions to its rate with `scale_time()` (`fn branches_at_rate()`, same file), so a probe at a higher rate has narrower branch-time distributions with a finer spacing than the run that just passed.

These backward passes run outside `pub(crate) fn run_timetree()` ([packages/treetime/src/timetree/inference/runner.rs#L44-L56](../../packages/treetime/src/timetree/inference/runner.rs#L44-L56)), which is the only caller of `fn explain_grid_point_limit()` ([packages/treetime/src/timetree/inference/runner.rs#L163-L204](../../packages/treetime/src/timetree/inference/runner.rs#L163-L204)). A grid point limit error of a probe therefore reaches the user without the explanation that [kb/decisions/distribution-grid-point-limit.md](../decisions/distribution-grid-point-limit.md) specifies: the `--max-grid-points` flag, the node name, the operands that set the span and the spacing, and the clock-rate advice.

## User-visible message

The message below follows from the wrappers on the error path, with placeholders in angle brackets; no run has produced it. The backward pass wraps the error with its grid request (`wrap_grid_request()` in [packages/treetime/src/timetree/inference/backward_pass.rs#L136](../../packages/treetime/src/timetree/inference/backward_pass.rs#L136), [#L162](../../packages/treetime/src/timetree/inference/backward_pass.rs#L162), [#L187](../../packages/treetime/src/timetree/inference/backward_pass.rs#L187)), whose display names the node by its graph key (`fn key_label()`, [packages/treetime/src/timetree/inference/grid_request.rs#L322-L324](../../packages/treetime/src/timetree/inference/grid_request.rs#L322-L324)) and not by its name. The probe, the round, and the refinement loop ([packages/treetime/src/timetree/refinement_loop.rs#L47-L48](../../packages/treetime/src/timetree/refinement_loop.rs#L47-L48)) add their contexts:

```text
Error:
   0: When running round <i>
   1: Failed to estimate the clock rate from the marginal likelihood of the node times
   2: When computing the marginal likelihood of the node times at clock rate <probe rate>
   3: When sending the time message of node <graph key> backward to its parent
   4: A grid over [<start>, <end>] with spacing <dx> needs <n> points, more than the limit of <limit>
```

For comparison, the same limit hit inside `run_timetree()` on `data/rsv/a/20` (`--coalescent opt --max-iter 2 --max-grid-points 8300`) gives:

```text
Error:
   0: Final timetree inference failed
   1: Time inference needs a grid of 8302 points, more than the limit of 8300 (--max-grid-points) when sending the time message forward to node 'NODE_0000009'. The grid spans [1931.44, 1991.11] at spacing 0.00718828: the span adds up the spans of the operands, and the spacing is that of the time distribution of the parent of node 'NODE_0000009' without the time message from node 'NODE_0000009', which takes it from the time distribution of the branch above node 'NODE_0000000'. The operands are [...]. The estimated clock rate 0.00061 may be far from the true rate, [...] Check the input dates, set --clock-rate, or raise --max-grid-points when enough memory is available.
   2: When sending the time message forward to node 16
   3: A grid over [1931.44, 1991.11] with spacing 0.00718828 needs 8302 points, more than the limit of 8300
```

The probe message lacks line 1 of this message. It names neither the flag nor the node name, and it does not say that the clock rate probe, not the inference at the current rate, needed the grid.

> [!IMPORTANT]
> **Investigation required.** No run has reproduced a probe error. 182 runs of `timetree --coalescent opt --max-iter 2` or `3` with `--max-grid-points` from 1000 up to the first passing value, on `rsv/a/20`, `flu/h3n2/200`, `ebola/100`, `ebola/362`, `zika/86`, `tb/149`, `lassa/L/200`, `dengue/100`, and `mpox/clade-ii/100`, all stopped in `run_timetree()` first, because the forward pass at the current rate needs larger grids than the backward pass of a probe. Find a run where a probe at a rising rate exceeds the limit while the inference at the current rate does not, for example a dataset whose rate step is large, and record the message.

## Fix direction

Explain the errors of the probes as `run_timetree()` explains its own, with the probe rate as the clock rate of the advice. The probe has no node names in `MarginalRateInputs`, so the explanation needs `NodeNames` passed in or applied by the caller in `fn marginal_clock_rate()`. Add a test that runs `step_marginal_clock_rate()` with a limit that only a probe exceeds and asserts the full message.
