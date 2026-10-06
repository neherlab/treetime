# Reading with several Newick dialects parses the whole tree once per rejected dialect

`fn read_tree_text()` in [packages/util-newick/src/read/select.rs](../../packages/util-newick/src/read/select.rs) tries the selected dialects in order. A dialect that defines `[&` comments reads a comment that does not follow its annotation syntax as a `malformed_annotation` token, in strict and in tolerant form alike ([packages/util-newick/src/newick.pest](../../packages/util-newick/src/newick.pest)). The grammar therefore accepts the whole tree, and only the mapping step rejects it at the first malformed annotation. Each rejected dialect costs one full parse.

Measured with `benches/read_dialects.rs` on a 2 MB annotated BEAST ladder tree of 20,000 leaves (bench profile, in the development container): 179 ms with `[Beast]`, 432 ms with `NewickDialect::ALL`, which first rejects Rich Newick, eNewick, NHX and MrBayes. TreeTime reads with `[Beast, Nhx, MrBayes, Classic]`, so BEAST and plain trees read in one parse; NHX trees cost one extra parse.

## Fix direction

Give the strict start rules comment rules without `malformed_annotation`, so a strict parse stops at the first `[&` comment that the dialect cannot read. pest rules take no parameters, so this duplicates the composite rules of the five dialects that reserve `[&` (comment choice, open, field, tail, items, preamble), and the strict and tolerant forms would differ in more than the start rule.

> [!IMPORTANT]
> **Decision required.** The duplicated rules make the grammar larger to save the cost of rejected dialects. Options: duplicate the strict composite rules, or keep one set of composite rules and accept one parse per rejected dialect. Evidence to weigh: the time of `NewickDialect::ALL` against a single dialect on large annotated trees of the expected inputs.
