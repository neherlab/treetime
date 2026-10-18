# Root committed at its median may be unneeded under the mean-field prior

Phase 3 of [timetree-coalescent-prior-and-clock-rate-streamlining.md](../proposals/timetree-coalescent-prior-and-clock-rate-streamlining.md).

The root is committed at its posterior median ([timetree-root-time-at-posterior-median.md](../decisions/timetree-root-time-at-posterior-median.md))
because its peak flipped between two modes on `data/ebola/20` under the old prior. Re-test under the mean-field prior:
return to the peak if the root is stable over rounds, otherwise keep the median and record why.
