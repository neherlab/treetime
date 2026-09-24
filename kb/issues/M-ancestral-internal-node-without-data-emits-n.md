# Internal nodes emit `N` where no descendant has data; v0 emits a base

At an alignment column where every descendant leaf of an internal node is `N`, v1 writes `N` for that internal node in both dense and sparse reconstruction. v0 writes the most likely base, which with no data is the most likely state under the model's equilibrium frequencies.

## Evidence

Tree `((L1,L2)X,L3)root` with columns 3-4 of every leaf set to `N`:

| Engine            | X    | root |
| ----------------- | ---- | ---- |
| v0 (`--gtr JC69`) | `CC` | `CC` |
| v1 dense          | `NN` | `NN` |
| v1 sparse         | `NN` | `NN` |

- v0 uses the 5-state `nuc` alphabet in which `N` is the all-ones profile ([packages/legacy/treetime/treetime/seq_utils.py#L20-L35](../../packages/legacy/treetime/treetime/seq_utils.py#L20-L35)) and assigns internal sequences by `prof2seq` argmax ([packages/legacy/treetime/treetime/treeanc.py#L909-L919](../../packages/legacy/treetime/treetime/treeanc.py#L909-L919)), so it never emits `N` at an internal node
- v1 marks such columns as unknown from the children (`fn compute_node_ranges` in [packages/treetime/src/seq/indel.rs](../../packages/treetime/src/seq/indel.rs)) and fills them with `N`

No decision records the divergence.

## Open question

Keep `N` (no data, so no inferred state) and record the decision, or emit the model's most likely base as v0 does.
