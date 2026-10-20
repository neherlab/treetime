# Polytomy resolution can stretch a time grid back to year 0 and stop the run

## Symptom

`treetime timetree` on `data/zika/86` with `--coalescent-opt --resolve-polytomies --seed 7 --max-iter 6` stops in
round 2:

```
When running round 2
Timetree inference failed
Time inference needs a grid of 1062601 points, more than the limit of 1000000 (--max-grid-points). The estimated clock rate 0.00111 may be far from the true rate, ...
A grid over [0, 2949.41] with spacing 0.00277566 needs 1062601 points, more than the limit of 1000000
```

The rate (1.1e-3) is plausible; a distribution reaching back to calendar year 0 is not, so the advice of the error
message about the clock rate does not apply here. The same command without `--resolve-polytomies` runs through and
writes all outputs.

## Investigation required

> [!IMPORTANT]
> **Investigation required.** Find the node or branch whose distribution spans to year 0 after the resolutions: a
> resolved node with a sampled time far in the past, a branch whose time length came out non-positive, or a
> polytomy-cleanup removal (see
> [M-timetree-polytomy-cleanup-drops-dated-single-child-node.md](M-timetree-polytomy-cleanup-drops-dated-single-child-node.md)).

## Validation

The reproduction completes with the default grid limit.

## Related issues

Other runs stop or slow down on the time grid through mass-window resampling and convolution at the finest operand spacing ([kb/proposals/timetree-forward-pass-grid-support.md](../proposals/timetree-forward-pass-grid-support.md)). Whether this run fails through the same mechanism is not shown.

- [M-timetree-node-support-reaches-far-past-at-moderate-clock-rate.md](M-timetree-node-support-reaches-far-past-at-moderate-clock-rate.md)
- [M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md](M-timetree-narrow-branch-distribution-forces-millions-of-grid-points.md)
- [M-timetree-grid-spacing-inherits-finest-operand.md](M-timetree-grid-spacing-inherits-finest-operand.md)
