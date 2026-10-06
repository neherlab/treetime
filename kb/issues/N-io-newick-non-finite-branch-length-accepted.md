# Newick reader accepts infinite branch lengths that no writer can write back

> [!IMPORTANT]
> **Decision required.** The reader accepts a branch length that overflows `f64`, such as `:1e400`, as infinity, as v0 does. The Newick writer rejects non-finite lengths, because `inf` is not a Newick number, so such a tree fails at its first tree output instead of at input. Options: reject non-finite lengths in the reader with the line and column of the length (a divergence from v0 that needs a `kb/decisions/` entry), or keep the v0 behavior.

## Evidence

- `fn Builder::set_field_value()` in [packages/util-newick/src/read/builder.rs](../../packages/util-newick/src/read/builder.rs) stores the parsed `f64` without a finiteness check
- The `util-newick` writer returns "When writing the branch above node 1 ('A'): Newick cannot represent the number inf" (`test_tree_read_infinite_length_cannot_be_written` in [packages/treetime-io/src/__tests__/test_tree_read.rs](../../packages/treetime-io/src/__tests__/test_tree_read.rs))
- Biopython 1.88, the v0 reader, reads `(A:1e400,B:1);` with branch length `inf`
