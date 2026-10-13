"""Run treetime on the simulated data. Resumable: a job with a successful run.json is skipped.

Usage: python run.py --suite timetree --set core -j 16   (suites: ancestral, timetree, optimize)
Results: results/<suite>/<variant>/<coal>/rep<NN>/L<L>_mu<mu>/
"""

import argparse
import json
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

from inputs import divergence_tree

HERE = Path(__file__).parent
COMMON = ['--no-progress', '-q', '-j', '1']
ENV = os.environ | {'OPENBLAS_NUM_THREADS': '1'}


def available_cpus():
    """CPUs this process may use: the cgroup quota when one is set, which os.cpu_count() does not see."""
    try:
        quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
        if quota != 'max':
            return max(1, int(quota) // int(period))
    except (OSError, ValueError):
        pass
    return len(os.sched_getaffinity(0))

# variant -> (input tree, extra flags)
ANCESTRAL = {
    'marginal': ('exact', []),
    'marginal-gtr3': ('exact', ['--gtr-iterations', '3']),
    'parsimony': ('exact', ['--method-anc', 'parsimony']),
}
KEEP_TOPOLOGY = ['--no-collapse-short-branches', '--no-merge-siblings', '--no-flip-parent-child']
OPTIMIZE = {
    'fixed-topo': ('perturbed', KEEP_TOPOLOGY),
    'fixed-topo-dense': ('perturbed', [*KEEP_TOPOLOGY, '--dense', 'true']),
    'default': ('perturbed', []),
}
ROOTING = {'keep': ('realized', ['--keep-root']), 'reroot': ('realized', ['--reroot', 'least-squares'])}
COALESCENT = {'none': [], 'opt': ['--coalescent-opt'], 'skyline': ['--coalescent-skyline']}
CLOCK = {'est': [], 'fixed': ['--clock-rate', '{mu}']}  # {mu} is the true rate of the sequence set
TIMETREE_CORE = {
    f'{r}-{c}-{k}': (tree, rflags + cflags + kflags)
    for r, (tree, rflags) in ROOTING.items()
    for c, cflags in COALESCENT.items()
    for k, kflags in CLOCK.items()
}
TIMETREE_EXTRA = {  # estimated clock rate
    'exact-keep-none-est': ('exact', ['--keep-root']),
    'rate-keep-none-est': ('rate', ['--keep-root']),
    'rate-reroot-opt-est': ('rate', ['--reroot', 'least-squares', '--coalescent-opt']),
    'keep-tc10-est': ('realized', ['--keep-root', '--coalescent', '10']),
    'keep-skyline-est-ci': ('realized', ['--keep-root', '--coalescent-skyline', '--time-marginal', 'only-final']),
}
EXTRA_SETS = r'/rep0[0-2]/L1000_mu0\.00[31]$'
SMOKE_SETS = r'/rep00/L1000_'


def command(suite, binary, tree, set_dir, flags, out):
    mu = json.loads((set_dir / 'gtr.json').read_text())['mu']
    flags = [f.format(mu=f'{mu:.10g}') for f in flags]
    aln = ['-t', str(tree), '-a', str(set_dir / 'tips.fasta.xz')]
    if suite == 'ancestral':
        outputs = ['--output-reconstructed-nuc-fasta', out / 'anc.fasta.xz']
        if 'parsimony' not in flags:  # parsimony fits no GTR model
            outputs += ['--output-gtr', out / 'gtr.json']
        return [binary, 'ancestral', *aln, *flags, *COMMON, *map(str, outputs)]
    if suite == 'optimize':
        outputs = ['--output-tree-nwk', out / 'tree.nwk', '--output-gtr', out / 'gtr.json']
        return [binary, 'optimize', *aln, *flags, *COMMON, *map(str, outputs)]
    dates = ['-d', str(set_dir.parent / 'tip_dates.tsv'), '--date-column', 'date', '--clock-filter', '0']
    outputs = ['--output-augur-node-data', out / 'node_data.json', '--output-tree-nwk', out / 'timetree.nwk',
               '--output-clock-model', out / 'clock.json']  # fmt: skip
    if any(f.startswith('--coalescent') for f in flags):  # requesting it without a model is an error
        outputs += ['--output-coalescent-tsv', out / 'coalescent.tsv']
    return [binary, 'timetree', *aln, *dates, *flags, *COMMON, *map(str, outputs)]


def run_job(job):
    cmd, out = job
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(out / 'log.txt', 'w') as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=ENV)
        _, status, usage = os.wait4(proc.pid, 0)
    rc = proc.returncode = os.waitstatus_to_exitcode(status)  # -9 is usually the OOM killer
    meta = {'cmd': cmd, 'returncode': rc, 'seconds': round(time.time() - t0, 1), 'max_rss_mb': usage.ru_maxrss // 1024}
    (out / 'run.json').write_text(json.dumps(meta, indent=1))
    print(f'{"ok " if rc == 0 else "FAIL"} {meta["seconds"]:7.1f}s {out}', flush=True)
    return rc


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--suite', choices=['ancestral', 'timetree', 'optimize'], required=True)
    p.add_argument('--set', choices=['smoke', 'core', 'extra'], default='smoke')
    p.add_argument('--filter', default='', help='regex on the job output path')
    p.add_argument('--binary', default=str(HERE.parent / 'target' / 'release' / 'treetime'))
    p.add_argument('--data', type=Path, default=HERE / 'data')
    p.add_argument('--results', type=Path, default=HERE / 'results')
    p.add_argument('-j', '--jobs', type=int, default=available_cpus())
    p.add_argument('--shard', default='0/1', help='I/N: run only every N-th sequence set, starting at I (Slurm arrays)')
    args = p.parse_args()
    shard, n_shards = map(int, args.shard.split('/'))

    variants = {'ancestral': ANCESTRAL, 'optimize': OPTIMIZE}.get(args.suite) or (
        TIMETREE_EXTRA if args.set == 'extra' else TIMETREE_CORE
    )
    sets = sorted(d for d in args.data.glob('*/rep*/L*') if (d / 'gtr.json').exists())
    sets = [d for d in sets if re.search({'smoke': SMOKE_SETS, 'core': '', 'extra': EXTRA_SETS}[args.set], str(d))]
    sets = sets[shard::n_shards]  # whole sets per shard: jobs of one set share its input trees

    plan = []
    for d in sets:
        for name, (tree, flags) in variants.items():
            out = args.results / args.suite / name / d.relative_to(args.data)
            if re.search(args.filter, str(out)):
                plan.append((d, tree, flags, out))

    scored = {'optimize': ['realized', 'exact']}.get(args.suite, [])  # truth trees that evaluate.py reads
    inputs = sorted({(d, tree) for d, tree, _, _ in plan} | {(d, kind) for d, *_ in plan for kind in scored})
    with ProcessPoolExecutor(args.jobs) as pool:  # input trees are shared between jobs, so write them first
        list(pool.map(divergence_tree, [d for d, _ in inputs], [t for _, t in inputs]))

    jobs = []
    for d, tree, flags, out in plan:
        done = out / 'run.json'
        if done.exists() and json.loads(done.read_text())['returncode'] == 0:
            continue
        jobs.append((command(args.suite, args.binary, divergence_tree(d, tree), d, flags, out), out))
    print(f'{len(jobs)} of {len(plan)} jobs to run', flush=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        failed = sum(rc != 0 for rc in pool.map(run_job, jobs))
    print(f'{failed} failed')


if __name__ == '__main__':
    main()
