# Type and reproduce performance reports

Build one typed benchmark-report pipeline that preserves every comparison dimension and reproduces all committed charts.

## Required changes

- Parse samples into a typed schema keyed by revision, workload, dataset, workers, and replicate.
- Reject missing or duplicate dimensions.
- Aggregate only after grouping by revision and the complete workload key.
- Benchmark the actual requested Rayon worker count.
- Generate paired, scaling, and CPU charts with the workspace Plotters stack.
- Regenerate every committed report from checked-in inputs.

## Validation

- Two-revision fixture that fails if revisions are averaged together.
- Worker-count fixture proving scaling samples use the requested pool.
- Golden structural tests for Markdown tables and chart series.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-benchmark-reports-mix-revisions-and-are-not-reproducible.md](../issues/M-benchmark-reports-mix-revisions-and-are-not-reproducible.md)
