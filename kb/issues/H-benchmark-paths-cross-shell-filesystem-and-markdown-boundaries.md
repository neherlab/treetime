# Benchmark paths cross shell, filesystem, and Markdown boundaries unsafely

Benchmark output and dataset strings are reused as shell fragments, filesystem paths, and generated Markdown without parsing for those domains.

The CLI validates paths by substring [dev/bench-graph-pass-cli#L86-L110](../../dev/bench-graph-pass-cli#L86-L110), interpolates the output path into Hyperfine's shell-valued `--prepare` command [dev/bench-graph-pass-cli#L232-L250](../../dev/bench-graph-pass-cli#L232-L250), and passes dataset labels into the report generator [dev/bench-graph-pass-report#L49-L54](../../dev/bench-graph-pass-report#L49-L54).

## Impact

- `--output` is interpolated into Hyperfine’s shell-valued `--prepare`, permitting command injection.
- Absolute or traversal-bearing dataset labels can write outside the selected output directory.
- Metadata can inject headings, links, or raw markup into generated reports.

## Potential solutions

- O1. Parse each boundary into a domain type and eliminate shell interpolation.
- O2. Escape raw strings independently for shell, path, and Markdown contexts. This retains three security-sensitive encoders and does not prove path containment.

## Recommendation

Eliminate the shell-valued cleanup command in favor of direct filesystem operations in the trusted harness. Parse output roots and dataset identifiers into bounded path types, verify normalized containment, and escape Markdown text at rendering time.

## Related issues

- [M-benchmark-reports-mix-revisions-and-are-not-reproducible.md](M-benchmark-reports-mix-revisions-and-are-not-reproducible.md)
