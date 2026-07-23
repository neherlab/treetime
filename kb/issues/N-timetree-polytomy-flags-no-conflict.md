# --keep-polytomies and --resolve-polytomies no conflicts_with declaration

`--keep-polytomies` and `--resolve-polytomies` are mutually exclusive flags but have no `conflicts_with` declaration in clap. Both can be passed simultaneously without error.

## Invalid state

The flags express competing topology policies: retain multifurcations or resolve them. Independent booleans permit both policies at once and force downstream code to apply implicit precedence. The same invalid representation can cross non-CLI boundaries when the options are copied into transport or pipeline parameters.

## Required contract

Represent polytomy handling as one exhaustive policy and parse CLI flags into it. If both legacy source flags are supplied, clap and every other transport adapter must reject the request with the same explanation. The internal pipeline must never receive contradictory booleans.

## Validation

- CLI parsing rejects both flags together and accepts each independently.
- Server and direct request parsing enforce the same state space.
- Existing keep and resolve outputs remain unchanged.
- An exhaustive match prevents a new policy variant from being silently ignored.

## Related tickets

- [kb/tickets/timetree-polytomy-flags-missing-conflicts-with-declaration.md](../tickets/timetree-polytomy-flags-missing-conflicts-with-declaration.md)

## Related issues

- [M-timetree-coalescent-modes-represent-invalid-states.md](M-timetree-coalescent-modes-represent-invalid-states.md)
