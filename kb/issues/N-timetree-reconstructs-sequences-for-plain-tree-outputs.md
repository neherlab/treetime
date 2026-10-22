# Timetree reconstructs sequences for outputs that do not use them

`sequence_outputs_requested()` ([packages/app-commands/src/commands/timetree/run.rs#L456-L462](../../packages/app-commands/src/commands/timetree/run.rs#L456-L462)) enables the final sequence reconstruction for every tree output except graph JSON and DOT. A Newick tree without `--newick-annotations` writes no mutations, so a run that asks only for plain Newick pays for a reconstruction it discards. Nexus outputs always carry the mutations in their BEAST comments, so they need the reconstruction.

Required behavior: enable the reconstruction only for outputs that read mutations or sequences: Auspice JSON, UShER MAT, Nexus, Newick with `--newick-annotations beast` or `nhx`, the reconstructed FASTA, and `--divergence-units=mutations`.
