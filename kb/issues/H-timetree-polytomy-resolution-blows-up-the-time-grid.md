# Polytomy resolution can stretch a time grid back to year 0 and stop the run

## Symptom

`treetime timetree` on `data/zika/86` with `--coalescent-opt --resolve-polytomies --seed 7 --max-iter 6` stops in
round 4 after resolving 5, 2, 2 and 2-3 polytomies:

```
Time inference needs a grid of 4141330 points, more than the limit of 1000000 (--max-grid-points).
A grid over [0, 9120.38] with spacing 0.00220228 needs 4141330 points, more than the limit of 1000000
```

The rate (1.4e-3) is plausible; a distribution reaching back to calendar year 0 is not. Both the phase-2 build and
the build whose prior follows topology changes fail the same way, so the cause predates both. Without
`--resolve-polytomies` the dataset runs through.

## Investigation required

> [!IMPORTANT]
> **Investigation required.** Find the node or branch whose distribution spans to year 0 after the resolutions: a
> resolved node with a sampled time far in the past, a branch whose time length came out non-positive, or a
> polytomy-cleanup removal (see
> [M-timetree-polytomy-cleanup-drops-dated-single-child-node.md](M-timetree-polytomy-cleanup-drops-dated-single-child-node.md)).

## Validation

The reproduction completes with the default grid limit.
