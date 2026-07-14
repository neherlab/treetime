# --keep-polytomies and --resolve-polytomies no conflicts_with declaration

`--keep-polytomies` and `--resolve-polytomies` are mutually exclusive flags but have no `conflicts_with` declaration in clap. Both can be passed simultaneously without error.

## Related tickets

- [kb/tickets/timetree-polytomy-flags-missing-conflicts-with-declaration.md](../tickets/timetree-polytomy-flags-missing-conflicts-with-declaration.md)
