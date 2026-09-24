# Dense ancestral and optimize runs need up to twice the peak memory of the rust branch

The dense (`--dense=true`) `ancestral` and `optimize` commands of the current development branch reach a peak resident memory 1.3 to 2.2 times that of the `rust` branch on the same inputs. Both binaries use jemalloc on Linux, so the growth comes from code changes, not from the allocator.

## Measurements

Peak resident set size (`ru_maxrss`) per process, measured by `dev/smoke` inside the build container. The `optimize` cases use `--dense=true`; the `ancestral` cases use `--method-anc=marginal --dense=true --model=jc69`.

| Case                                          |   `rust` |          current |
| --------------------------------------------- | -------: | ---------------: |
| `optimize/mpox/clade-ii/100/dense`            |  5.3 GiB |         11.0 GiB |
| `optimize/mpox/clade-ii/500/dense`            | 22.8 GiB |         49.9 GiB |
| `ancestral/mpox/clade-ii/500/marginal-dense`  | 22.6 GiB |         33.6 GiB |
| `optimize/mpox/clade-ii/1000/dense`           | 46.9 GiB | more than 78 GiB |
| `ancestral/mpox/clade-ii/1000/marginal-dense` | 46.7 GiB | more than 78 GiB |
| `optimize/sc2/4500/dense`                     | 38.1 GiB | more than 78 GiB |

"More than 78 GiB" means the process failed with an allocation error at the `dev/smoke` memory limit of 78 GiB, which the `rust` binary stayed below.

## Impact and scope

- Dense runs on the larger `mpox/clade-ii` and `sc2` datasets no longer fit on a 128 GB machine next to other work.
- `ancestral/mpox/clade-ii/2000/marginal-dense` and `optimize/mpox/clade-ii/2000/dense` exceed the limit on both branches.

## Reproduction

```bash
./dev/smoke --tier full --only 'mpox/clade-ii/(100|500|1000)/(dense|marginal-dense)$|sc2/4500/dense$'
```

The peak of each run is `max_rss_kb` in `snapshots/<id>/cases/<case-id>/result.json`.

## Related issues

- [N-io-large-dataset-memory-constraint.md](N-io-large-dataset-memory-constraint.md)
- [N-ancestral-marginal-array-kernels-allocate.md](N-ancestral-marginal-array-kernels-allocate.md)
- [M-timetree-marginal-dense-mpox-slow.md](M-timetree-marginal-dense-mpox-slow.md)
