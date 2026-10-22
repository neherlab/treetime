"""Profile of TreeTime's marginal log-likelihood of node times over fixed clock rates, from separate CLI runs.

Purpose
  For one `treetime` binary and one simulated sequence set of `validation/`, run `treetime timetree` once per clock
  rate of a geometric grid, with and without a Kingman coalescent prior, and collect `log_lh_time` of the trace log.
  The maximum of the profile over the rate shows where a run at a fixed rate fits the data best; profiles of two
  binaries on the same set show how a code change moves it.

Model
  Each run uses TreeTime's own model: branch likelihoods from the alignment, node times integrated on grids, and,
  with `--coalescent <Tc>`, the mean-field coalescent prior, whose expected lineage counts are frozen from the node
  posteriors of the same run. `log_lh_time` is the log marginal likelihood of the node times that the run reports
  after its single round (`--max-iter 1`) at the given `--clock-rate`.

  Each point builds its own prior at its own rate. The profile is therefore the end-to-end result of separate runs,
  and it differs from the objective of TreeTime's rate step, which rescales the branch likelihoods of one run and
  keeps that run's frozen prior. It differs as well from `poisson_timetree_surface.py`, which uses Poisson branch
  likelihoods from realized substitution counts and the exact Kingman density of the coupled node times.

Inputs (layout of `validation/README.md`)
  --data-set validation/data/<coal>/repNN/L<L>_mu<mu>
    tree       <data-set>/inputs/realized.nwk (true topology, branch lengths = realized differences / L). When it is
               missing, the script writes it with `divergence_tree()` of `validation/inputs.py`, as `validation/run.py`
               does; this step needs numpy and Biopython in `--python`
    alignment  <data-set>/tips.fasta.xz
    dates      <data-set>/../tip_dates.tsv, column `date`
  The command has the shape of `command()` in `validation/run.py` (timetree suite, `COMMON` flags):
    treetime timetree -t <tree> -a <aln> -d <dates> --date-column date --clock-filter 0 --clock-rate <mu> \
      --max-iter 1 --keep-root [--coalescent <Tc>] --no-progress -q -j 1 --output-tracelog <csv>
  with OPENBLAS_NUM_THREADS=1.

Output
  TSV with one row per run: data_set, true_mu (from gtr.json), mu, prior (`coalescent=<Tc>` or `none`),
  log_lh_time, log_lh_coal and log_lh_seq of the last trace-log row, rounds (rows of the trace log), seconds, error
  (empty on success), treetime_version (`treetime --version`), binary_sha256. A failed run keeps its row with the last
  line of its error output; the script then exits with status 1 after it writes the table.

Options
  --binary PATH      treetime executable (default: .out/treetime)
  --data-set DIR     sequence set directory (required)
  --tree PATH        input tree instead of <data-set>/inputs/realized.nwk
  --rates MIN MAX N  geometric rate grid (default: 5e-5 1.5e-3 15)
  --tc TC            coalescent time scale of the prior runs (default: 10)
  --output PATH      TSV file (required)
  --jobs N           runs in parallel (default: 1); each run uses one thread and up to about 2 GB at 400 tips
  --python PATH      Python that writes a missing input tree (default: this interpreter)

Examples (from the repository root; the binary needs the C library of the build image or the host)
  python test_scripts/timetree_log_z_profile.py --binary .out/treetime \
    --data-set validation/data/const/rep00/L1000_mu0.0003 --jobs 10 --output tmp/log_z_profile/rep00.tsv
  python test_scripts/timetree_log_z_profile.py --data-set validation/data/const/rep03/L10000_mu0.001 \
    --rates 2e-4 4e-3 21 --output tmp/log_z_profile/rep03_L10000.tsv

The script uses only the Python standard library, so it also runs with the system Python of the build image.
"""

import argparse
import csv
import hashlib
import json
import logging
import math
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
VALIDATION = REPO / 'validation'
COMMON = ['--no-progress', '-q', '-j', '1']
ENV = os.environ | {'OPENBLAS_NUM_THREADS': '1'}
TRACE_COLUMNS = ('log_lh_time', 'log_lh_coal', 'log_lh_seq')

log = logging.getLogger('timetree_log_z_profile')


def main() -> None:
  logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stderr)
  args = parse_args()
  try:
    data = DataSet.of(args.data_set, args.tree, args.python)
    binary = BinaryIdentity.of(args.binary)
  except (OSError, subprocess.CalledProcessError) as error:
    log.error(f'error: {error}')
    sys.exit(1)
  rates = geometric_grid(*args.rates)
  priors = [Prior(args.tc), Prior(None)]
  runs = [(mu, prior) for prior in priors for mu in rates]
  log.info(f'{binary.version} ({binary.sha256[:12]}), {data.set_dir}: {len(runs)} runs on {args.jobs} threads')
  with tempfile.TemporaryDirectory(prefix='log_z_profile_') as work, ThreadPoolExecutor(args.jobs) as pool:
    rows = list(pool.map(lambda job: run_point(args.binary, data, job[0], job[1], Path(work), binary), runs))
  write_table(args.output, rows)
  failed = [row for row in rows if row.error]
  if failed:
    log.error(f'{len(failed)} of {len(rows)} runs failed; see the error column of {args.output}')
    sys.exit(1)


def parse_args() -> argparse.Namespace:
  p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  p.add_argument('--binary', type=Path, default=REPO / '.out' / 'treetime')
  p.add_argument('--data-set', type=Path, required=True)
  p.add_argument('--tree', type=Path)
  p.add_argument('--rates', type=float, nargs=3, default=[5e-5, 1.5e-3, 15], metavar=('MIN', 'MAX', 'N'))
  p.add_argument('--tc', type=float, default=10.0)
  p.add_argument('--output', type=Path, required=True)
  p.add_argument('--jobs', type=int, default=1)
  p.add_argument('--python', type=Path, default=Path(sys.executable))
  args = p.parse_args()
  if not 0 < args.rates[0] < args.rates[1] or args.rates[2] < 2 or args.rates[2] != int(args.rates[2]):
    p.error('--rates needs 0 < MIN < MAX and an integer N >= 2')
  return args


def geometric_grid(low: float, high: float, n: float) -> list[float]:
  count = int(n)
  return [low * (high / low) ** (i / (count - 1)) for i in range(count)]


@dataclass(frozen=True)
class BinaryIdentity:
  version: str
  sha256: str

  @classmethod
  def of(cls, binary: Path) -> 'BinaryIdentity':
    if not binary.is_file():
      raise FileNotFoundError(f'treetime binary not found: {binary}')
    version = subprocess.run([binary, '--version'], capture_output=True, text=True, check=True).stdout.strip()
    digest = hashlib.sha256()
    with binary.open('rb') as fh:
      for chunk in iter(lambda: fh.read(1 << 20), b''):
        digest.update(chunk)
    return cls(version, digest.hexdigest())


@dataclass(frozen=True)
class DataSet:
  set_dir: Path
  tree: Path
  alignment: Path
  dates: Path
  true_mu: float

  @classmethod
  def of(cls, set_dir: Path, tree: Path | None, python: Path) -> 'DataSet':
    gtr = set_dir / 'gtr.json'
    if not gtr.is_file():
      raise FileNotFoundError(f'not a sequence set of validation/ (no gtr.json): {set_dir}')
    true_mu = json.loads(gtr.read_text())['mu']
    return cls(set_dir, tree or realized_tree(set_dir, python), set_dir / 'tips.fasta.xz',
               set_dir.parent / 'tip_dates.tsv', true_mu)


def realized_tree(set_dir: Path, python: Path) -> Path:
  path = set_dir / 'inputs' / 'realized.nwk'
  if path.exists():
    return path
  log.info(f'writing {path} with divergence_tree() of validation/inputs.py')
  script = 'import sys; from pathlib import Path; from inputs import divergence_tree; ' \
           'divergence_tree(Path(sys.argv[1]), "realized")'
  subprocess.run([python, '-c', script, set_dir.resolve()], cwd=VALIDATION, check=True)
  return path


@dataclass(frozen=True)
class Prior:
  tc: float | None

  def flags(self) -> list[str]:
    return ['--coalescent', f'{self.tc:g}'] if self.tc is not None else []

  def label(self) -> str:
    return f'coalescent={self.tc:g}' if self.tc is not None else 'none'


def run_point(binary: Path, data: DataSet, mu: float, prior: Prior, work: Path,
              identity: BinaryIdentity) -> 'ProfileRow':
  trace = work / f'trace_{prior.label()}_{mu:.10g}.csv'
  command = [
    binary, 'timetree', '-t', data.tree, '-a', data.alignment,
    '-d', data.dates, '--date-column', 'date', '--clock-filter', '0',
    '--clock-rate', f'{mu:.10g}', '--max-iter', '1', '--keep-root', *prior.flags(),
    *COMMON, '--output-tracelog', trace,
  ]  # fmt: skip
  started = time.time()
  proc = subprocess.run(list(map(str, command)), capture_output=True, text=True, env=ENV)
  seconds = round(time.time() - started, 1)
  values, rounds, error = read_trace(trace) if proc.returncode == 0 else ({}, 0, last_line(proc))
  log.info(f'{"ok  " if not error else "FAIL"} {seconds:6.1f}s mu={mu:.4e} {prior.label()} '
           f'log_lh_time={values.get("log_lh_time", math.nan):.4f}')
  return ProfileRow(
    data_set=str(data.set_dir), true_mu=data.true_mu, mu=mu, prior=prior.label(),
    log_lh_time=values.get('log_lh_time', math.nan), log_lh_coal=values.get('log_lh_coal', math.nan),
    log_lh_seq=values.get('log_lh_seq', math.nan), rounds=rounds, seconds=seconds, error=error,
    treetime_version=identity.version, binary_sha256=identity.sha256,
  )


def read_trace(trace: Path) -> tuple[dict[str, float], int, str]:
  if not trace.is_file():
    return {}, 0, 'the run wrote no trace log'
  with trace.open(newline='') as fh:
    records = list(csv.DictReader(fh))
  if not records:
    return {}, 0, 'the trace log has no rows'
  return {column: float(records[-1][column] or math.nan) for column in TRACE_COLUMNS}, len(records), ''


def last_line(proc: subprocess.CompletedProcess) -> str:
  lines = [line for line in (proc.stderr + proc.stdout).splitlines() if line.strip()]
  return f'exit {proc.returncode}: {lines[-1].strip() if lines else "no output"}'


@dataclass(frozen=True)
class ProfileRow:
  data_set: str
  true_mu: float
  mu: float
  prior: str
  log_lh_time: float
  log_lh_coal: float
  log_lh_seq: float
  rounds: int
  seconds: float
  error: str
  treetime_version: str
  binary_sha256: str


def write_table(path: Path, rows: list[ProfileRow]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open('w', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=list(ProfileRow.__dataclass_fields__), delimiter='\t',
                            lineterminator='\n')
    writer.writeheader()
    writer.writerows({k: format_value(v) for k, v in asdict(row).items()} for row in rows)


def format_value(value: object) -> object:
  return f'{value:.10g}' if isinstance(value, float) else value


if __name__ == '__main__':
  main()
