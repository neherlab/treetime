# Node dates TSV output unimplemented

> [!WARNING]
> **Needs review.** The confidence-interval variant of v0's `dates.tsv` exists. `--output-confidence-tsv` writes `timetree.confidence.tsv` ("Date intervals of every node", [`output_plan.rs#L415`](../../packages/app-output/src/output_plan.rs#L415), [`output_plan.rs#L478`](../../packages/app-output/src/output_plan.rs#L478)) with the columns `name`, `date` (numeric date), `lower`, `upper` ([`confidence.rs#L134-L142`](../../packages/treetime/src/timetree/confidence.rs#L134-L142)). The pipeline computes these intervals only with `--time-marginal=only-final`, `--time-marginal=always`, or a rate standard deviation ([`pipeline.rs#L471-L474`](../../packages/treetime/src/timetree/pipeline.rs#L471-L474)). What remains missing is a node dates table for a default run, and v0's calendar-date column next to the numeric date.

No `dates.tsv` output file is produced for a default run. v0 writes `dates.tsv` with columns including node name, date estimate, and (with `--confidence --covariation`) lower/upper bounds.

The same node date data (node name, `numdate`, resolved `date`, and `num_date_confidence` lower/upper bounds) is emitted in `timetree.augur-node-data.json`.

The public tracker documents `dates.tsv` as the Python command's textual confidence-interval output, including node, date, numeric date, and 90% posterior-region bounds [[issue](https://github.com/neherlab/treetime/issues/64)] [[comment](https://github.com/neherlab/treetime/issues/64#issuecomment-416872300)].
