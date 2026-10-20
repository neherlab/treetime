# DOT output writes an invisible edge for every pair of leaves

`print_fake_edges()` writes an invisible edge for every ordered pair of roots and of leaves, to keep them in one row ([packages/treetime-io/src/graphviz.rs#L135-L154](../../packages/treetime-io/src/graphviz.rs#L135-L154), called at [#L67](../../packages/treetime-io/src/graphviz.rs#L67) and [#L79](../../packages/treetime-io/src/graphviz.rs#L79)). 2000 leaves give 4 million lines.

Required behavior: stop writing these edges. Delete `print_fake_edges()` and its weight constants.
