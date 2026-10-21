# Mugration probabilities are written with three precisions

One per-node state profile is written three ways: the confidence CSV writes every state with 6 decimals ([packages/app-output/src/trait_tables.rs#L46](../../packages/app-output/src/trait_tables.rs#L46)), augur node data writes states above 0.001 at full precision, and Auspice JSON writes states above 0.001 rounded to 3 significant digits, and the entropy too ([packages/app-output/src/auspice.rs#L307-L316](../../packages/app-output/src/auspice.rs#L307-L316); the 0.001 filter is `build_confidence_map()`, [packages/app-output/src/trait_profile.rs#L9-L14](../../packages/app-output/src/trait_profile.rs#L9-L14)). All three read the same `StatePosterior` of the node.

Required behavior: every output writes probabilities and the entropy at full precision. The state filters stay as they are.
