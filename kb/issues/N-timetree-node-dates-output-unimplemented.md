# Node dates TSV output unimplemented

A default timetree run writes no node dates table. v0 writes `dates.tsv` on every timetree run, with the columns `node`, `date` (calendar date), and `numeric date`; with `--confidence` and either `--covariation` or `--clock-std-dev`, it adds the lower and upper bounds of the 90% max posterior region ([packages/legacy/treetime/treetime/CLI_io.py#L135-L157](../../packages/legacy/treetime/treetime/CLI_io.py#L135-L157), [packages/legacy/treetime/treetime/wrappers.py#L452-L467](../../packages/legacy/treetime/treetime/wrappers.py#L452-L467)).

## Current v1 output

The confidence-interval variant exists. `--output-confidence-tsv` writes `timetree.confidence.tsv` ("Date intervals of every node", [packages/app-output/src/output_plan.rs#L836](../../packages/app-output/src/output_plan.rs#L836)) with the columns `name`, `date` (numeric date), `lower`, `upper` ([packages/treetime/src/timetree/confidence.rs#L176-L184](../../packages/treetime/src/timetree/confidence.rs#L176-L184)). This output has two limits:

- The pipeline computes the intervals only with `--time-marginal=only-final`, `--time-marginal=always`, or a rate standard deviation ([packages/treetime/src/timetree/pipeline.rs#L453-L462](../../packages/treetime/src/timetree/pipeline.rs#L453-L462))
- `--output-all` without `--output-selection` leaves it out, because it is in the list of outputs that are not written by default ([packages/app-output/src/output_plan.rs#L605-L612](../../packages/app-output/src/output_plan.rs#L605-L612))

The same node date data (node name, `numdate`, calendar `date`, and `num_date_confidence` lower and upper bounds) is in `timetree.augur-node-data.json`.

## Remaining work

- A node dates table for a default run, with name and numeric date of every node
- v0's calendar-date column next to the numeric date

The public tracker documents `dates.tsv` as the Python command's textual confidence-interval output, including node, date, numeric date, and 90% posterior-region bounds [[issue](https://github.com/neherlab/treetime/issues/64)] [[comment](https://github.com/neherlab/treetime/issues/64#issuecomment-416872300)].
