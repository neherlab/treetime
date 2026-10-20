# Timetree reconstructs sequences for outputs that do not use them

`sequence_outputs_requested()` ([packages/app-commands/src/commands/timetree/run.rs#L470-L476](../../packages/app-commands/src/commands/timetree/run.rs#L470-L476)) enables the final sequence reconstruction for every tree output except graph JSON and DOT. Plain-style Newick and Nexus write no mutations, so a run that asks only for them pays for a reconstruction it discards.

Required behavior: enable the reconstruction only for outputs that read mutations or sequences: Auspice JSON, UShER MAT, Newick and Nexus in an annotated style (`--output-nwk-style=beast` or `nhx`), the reconstructed FASTA, and `--divergence-units=mutations`.
