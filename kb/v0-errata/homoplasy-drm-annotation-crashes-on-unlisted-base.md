# Homoplasy DRM annotation crashes on a derived base the table does not list

## v0 location

`read_in_DRMs()` keys the drug resistance table by position and stores the amino-acid substitution of each listed `ALT_BASE` [`packages/legacy/treetime/treetime/CLI_io.py#L44-L70`](../../packages/legacy/treetime/treetime/CLI_io.py#L44-L70). The ranked lists of `scan_homoplasies()` annotate a mutation at a listed position with `drms[mut[1]]['alt_base'][mut[2]]` [`packages/legacy/treetime/treetime/wrappers.py#L243`](../../packages/legacy/treetime/treetime/wrappers.py#L243) [`#L277`](../../packages/legacy/treetime/treetime/wrappers.py#L277).

## Erratum

The lookup checks the position (`if mut[1] in drms`) but not the derived base. A recurrent mutation at a listed position to a base that the table does not list raises `KeyError`, and the run stops after printing the header of the ranked list.

## Evidence

- The table format has one row per position and alternative base, so a position with some bases listed and others not is normal input. For example, a table that lists only `9344 T` for `data/zika/86` makes `treetime homoplasy --gtr jc --drms <table>` exit with `KeyError: np.str_('A')`, because `G9344A` occurs on four branches
- The adjacent per-taxon DRM count tests the position only (`mut[1] in drms`) [`packages/legacy/treetime/treetime/wrappers.py#L307`](../../packages/legacy/treetime/treetime/wrappers.py#L307), so the same mutation counts as a DRM there without a base lookup

## v0 impact

A DRM table that lists some alternative bases of a position stops every run with a recurrent mutation at that position to another base.

## v1 status

v1 annotates a substitution at a listed position with the gene and the drug, and adds the amino-acid substitution only when the table lists its derived base (`fn DrmTable::annotate()` in [`packages/app-commands/src/commands/homoplasy/drms.rs`](../../packages/app-commands/src/commands/homoplasy/drms.rs)). The table above gives `G9344A	4	NS5 TESTDRUG`, as recorded in [kb/decisions/homoplasy-mutation-mapping-and-counting.md](../decisions/homoplasy-mutation-mapping-and-counting.md).
