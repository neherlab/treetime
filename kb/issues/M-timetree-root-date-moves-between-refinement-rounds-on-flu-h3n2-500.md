# The root date of `flu/h3n2/500` moves by about 0.45 y between refinement rounds

On `flu/h3n2/500`, the root date of `treetime timetree` changes from round to round of the refinement loop, independent of the seed. After the default 2 rounds the root is at 1964.18 for every seed; with `--keep-polytomies` it is at 1964.64, and after 8-10 rounds at 1964.62. The shift of about 0.45 y is about 1000 times the spread over seeds (sd 0.00045 y). It lies inside the reported 90% root interval (1962.654 to 1965.653, 3.0 y), so the interval covers it.

## Evidence

Measured with `treetime timetree -v --seed <s>` at commit `1e8fbcd8b`, default settings, release build, seeds 1-20 and controls:

- `max_dt` (the largest node-time change of a round) stays at 0.87-1.03 y from round 2 to round 7, then falls below 0.01 y
- Rounds that resolve polytomies move the root far; `ebola/100` and `mpox/clade-ii/500` show no such move (root sd 0.0049 y and 0.0068 y over seeds)
- The move does not depend on the seed: every seed gives 1964.18 after 2 rounds

> [!IMPORTANT]
> **Investigation required.** The cause is unknown. Trace which node moves by about 1 y per round (node dates of seed 1 after 2 and after 10 rounds, against the run with `--keep-polytomies`), whether the polytomy resolution or the time inference of the next round moves it, and run v0 on the same dataset with `--max-iter 10` to check whether v0 moves the root the same way.

## Validation

- The cause is known, and either a fix makes the root date stable across rounds or a recorded explanation shows that the move is the expected effect of the resolution
