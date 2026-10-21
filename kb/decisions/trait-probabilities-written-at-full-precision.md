# Trait probabilities and entropy are written at full precision

Every output of a reconstructed trait writes the state probabilities and the entropy of a node at full precision, as the shortest decimal that reads back as the same double. The Auspice JSON therefore diverges from `augur export v2`, which writes them with 3 significant digits.

**Type**: Output of a scientific value; divergence from augur.

**Status**: Approved on 2026-10-21.

## Rule

- **Confidence CSV**: every state of every node, written through the table serializer like every other table ([packages/app-output/src/trait_tables.rs](../../packages/app-output/src/trait_tables.rs))
- **augur node data**: the states above 0.001 and the entropy ([packages/app-output/src/augur_node_data_traits.rs](../../packages/app-output/src/augur_node_data_traits.rs))
- **Auspice JSON**: the states above 0.001 and the entropy of each trait attribute ([packages/app-output/src/auspice.rs](../../packages/app-output/src/auspice.rs))
- The filter of states above 0.001 stays where augur applies it (`fn build_confidence_map()`, [packages/app-output/src/trait_annotations.rs](../../packages/app-output/src/trait_annotations.rs)). Date confidence intervals in the Auspice JSON keep the 3 significant digits of augur

All three outputs read the same `StatePosterior` of the node.

## Reference behavior

- **v0**, confidence CSV: every probability with `str(x)`, the shortest decimal that reads back ([packages/legacy/treetime/treetime/wrappers.py#L905-L911](../../packages/legacy/treetime/treetime/wrappers.py#L905-L911)). v1 matches
- **augur traits**, node data: the profile entries above 0.001 and the entropy, unrounded (`augur/traits.py` lines 84-95 at revision `58f329042e7fc631ca275cf2df44c8720f80949f`). v1 matches
- **augur export v2**, Auspice JSON: the confidences above 0.001 and the entropy with 3 significant digits, "to minimise the JSON size" (`format_number()` and `attr_confidence()`, `augur/export_v2.py` lines 689-699 and 742-750 at the same revision). v1 diverges

## Why

- One posterior gives one value in every output, so values read from the Auspice JSON, the CSV, and augur node data agree and can be compared without a rounding difference
- Three significant digits make close states look equal: posteriors 0.39165 and 0.39151 both show as 0.392, so a reader of the Auspice JSON cannot see which state the reconstruction chose
- The size cost is small: a node writes only the states above 0.001, so each node gains at most a few numbers of up to 17 significant digits
