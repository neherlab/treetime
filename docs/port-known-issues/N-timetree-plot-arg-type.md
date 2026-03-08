# --plot-rtt and --plot-tree typed as Option\<usize\>

`--plot-rtt` and `--plot-tree` are typed as `Option<usize>` but should be
`Option<PathBuf>` for output path specification.

## Related issues

- [Plot commands unimplemented](N-timetree-plot-unimplemented.md) the flags
  also panic with `todo!()`
