# Mugration probabilities are written with three precisions

One per-node state profile is written three ways: the confidence CSV writes every state with 6 decimals ([packages/app-output/src/trait_tables.rs#L46](../../packages/app-output/src/trait_tables.rs#L46)), augur node data writes states above 0.001 at full precision, and Auspice JSON writes states above 0.001 rounded to 3 significant digits, and the entropy too ([packages/app-output/src/auspice.rs#L307-L316](../../packages/app-output/src/auspice.rs#L307-L316); the 0.001 filter is `build_confidence_map()`, [packages/app-output/src/trait_annotations.rs#L9-L14](../../packages/app-output/src/trait_annotations.rs#L9-L14)). All three read the same `StatePosterior` of the node.

## Reference behavior

- **Confidence CSV**: v0 writes every probability with `str(x)`, the shortest decimal that reads back to the same float ([packages/legacy/treetime/treetime/wrappers.py#L905-L911](../../packages/legacy/treetime/treetime/wrappers.py#L905-L911)). The 6 decimals of v1 diverge from v0
- **augur node data**: `augur traits` stores the profile entries above 0.001 and the entropy unrounded (`augur/traits.py` lines 84-95 at revision `58f329042e7fc631ca275cf2df44c8720f80949f`). v1 matches
- **Auspice JSON**: `augur export v2` writes the confidences above 0.001 and the entropy with 3 significant digits, "to minimise the JSON size" (`format_number()` and `attr_confidence()`, `augur/export_v2.py` lines 689-699 and 742-750 at the same revision). v1 matches

> [!IMPORTANT]
> **Decision required.** Full precision in every output, as this issue asked, diverges from `augur export v2` in the Auspice JSON. Options:
>
> - Full precision in the confidence CSV only (v0 parity); augur node data and Auspice JSON stay as they are (augur parity)
> - Full precision in every output, and record the divergence of the Auspice JSON from augur
