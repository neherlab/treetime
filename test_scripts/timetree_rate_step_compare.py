"""Compare variants of the timetree rate step on the validation suite, paired with a baseline variant.

Purpose
  Summarize the timetree `core` runs with a coalescent prior and an estimated rate (`*-opt-est`, `*-skyline-est`) of
  several builds or settings of the rate step: rate error |ln(rate / true)|, node date RMSE, root date error, the
  share of sets where a variant has a smaller or larger rate error than the baseline, failures, and the share of runs
  within 10% and beyond 2-fold of the true rate.

Input
  A folder with one results folder per variant, each holding the `timetree.tsv` of `validation/evaluate.py`. The
  research build of branch `research/rate-step` reads the step settings from TREETIME_RATE_COORD (log, sqrt, cbrt,
  rate), TREETIME_RATE_GUARD (none: follow the parabola maximum in every round), TREETIME_RATE_STEP (step fraction)
  and TREETIME_RATE_BOUND (largest factor of one target); in validation/, per variant:
    pixi run python run.py --suite timetree --set core --binary <build> --results <root>/<variant> \
      --filter='-(opt|skyline)-est/'
    pixi run python evaluate.py --suite timetree --results <root>/<variant>

Usage
  pixi run python test_scripts/timetree_rate_step_compare.py <root> [baseline variant, default: base]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(sys.argv[1])
base = sys.argv[2] if len(sys.argv) > 2 else 'base'
KEY = ['variant', 'coal', 'rep', 'L', 'mu']
frames = []
for d in sorted(p for p in root.iterdir() if (p / 'timetree.tsv').exists()):
    t = pd.read_csv(d / 'timetree.tsv', sep='\t')
    frames.append(t.assign(step=d.name))
t = pd.concat(frames)
t['ok'] = t.returncode == 0
t['abs_log_rate'] = np.abs(np.log1p(t.clock_rate_rel_error))
t['abs_root'] = t.root_date_error.abs()
b = t[t.step == base].set_index(KEY)
t = t.join(b[['abs_log_rate', 'date_rmse']].add_suffix('_base'), on=KEY)
t['rate_better'] = t.abs_log_rate < t.abs_log_rate_base - 1e-9
t['rate_worse'] = t.abs_log_rate > t.abs_log_rate_base + 1e-9
pd.set_option('display.width', 250)
order = [base] + [s for s in t.step.unique() if s != base]
for metric, agg in [('abs_log_rate', 'median'), ('date_rmse', 'median'), ('abs_root', 'median'),
                    ('rate_better', 'mean'), ('rate_worse', 'mean')]:  # fmt: skip
    print(f'\n## {metric} ({agg}), all four variants pooled')
    p = t[t.ok].pivot_table(index=['L', 'mu'], columns='step', values=metric, aggfunc=agg)[order]
    p.loc[('all', ''), :] = t[t.ok].groupby('step')[metric].agg(agg)[order]
    print(p.round(3).to_string())
print('\n## failed runs')
print(t.assign(fail=~t.ok).pivot_table(index=['L', 'mu'], columns='step', values='fail', aggfunc='sum')[order].to_string())
print('\n## rate within 10% / beyond 2-fold of the truth (fraction of runs)')
for name, sel in [('within 10%', t.abs_log_rate < np.log(1.1)), ('beyond 2-fold', t.abs_log_rate > np.log(2))]:
    print(name)
    print(t.assign(s=sel).pivot_table(index=['L', 'mu'], columns='step', values='s', aggfunc='mean')[order].round(3).to_string())
print('\n## median abs log rate by timetree variant')
print(t[t.ok].pivot_table(index='variant', columns='step', values='abs_log_rate', aggfunc='median')[order].round(3).to_string())
