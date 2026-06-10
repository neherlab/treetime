# Newick writer emits per-key annotation blocks instead of standard BEAST format

`nwk_write_with` at `packages/treetime-io/src/nwk.rs#L254` writes each annotation as a separate `[&k="v"]` block:

```rust
.map(|(key, val)| format!("[&{key}=\"{val}\"]"))
```

Output: `name:0.1[&mutations="G42T"][&date="2020.5"]`

BEAST/FigTree standard: `name:0.1[&mutations="G42T",date="2020.5"]` -- all annotations in one block, comma-separated.

No major tool writes per-key blocks. Most tools would parse only the first block and ignore or reject subsequent ones. Biopython, ETE, DendroPy, FigTree, IQ-TREE all write a single `[&k1=v1,k2=v2]` block.

## Impact

Annotated trees written by treetime may not be correctly read by BEAST, FigTree, or other standard phylogenetic visualization tools.

## Locations

- `packages/treetime-io/src/nwk.rs#L250-L259`

## Related issues

- [kb/tickets/io-nwk-writer-style-dispatch.md](../tickets/io-nwk-writer-style-dispatch.md)
