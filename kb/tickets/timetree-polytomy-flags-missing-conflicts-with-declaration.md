# --keep-polytomies and --resolve-polytomies no conflicts_with declaration

`--keep-polytomies` and `--resolve-polytomies` are mutually exclusive flags but have no `conflicts_with` declaration in clap. Both can be passed simultaneously without error.

## Related issues

- Source: [kb/issues/N-timetree-polytomy-flags-no-conflict.md](../issues/N-timetree-polytomy-flags-no-conflict.md) -- delete after full resolution
