# Command outcome lists output files left by earlier runs

`CommandOutcome.output_files` tells the app which files a run wrote. `fn written_files()` [packages/app-commands/src/command.rs#L423-L431](../../packages/app-commands/src/command.rs#L423-L431) builds the list after the run from the output plan and the file system, not from the files the run wrote:

- A planned path is listed when a file exists there. An output that the run skips because it has no data (for example the GTR file under `--output-all` without a fitted model) is listed when an earlier run into the same folder wrote it
- A per-CDS amino-acid FASTA template is listed by a folder scan in `fn existing_files()` [packages/app-commands/src/command.rs#L433-L464](../../packages/app-commands/src/command.rs#L433-L464): every file whose name starts with the text before the first placeholder and ends with the text after it. Files of CDSes that the run did not reconstruct match too
- The scan splits the file name at the first placeholder only, so a template with both placeholders (`aa.{cds}.%GENE.fasta`) matches no written file

The run knows the exact paths: the output plan and, for the amino-acid FASTA, the expanded per-CDS paths that the ancestral runner checks for collisions (`fn cds_output_paths()` in [packages/app-commands/src/commands/ancestral/aa_node_data.rs](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs)).

## Fix direction

Return the written paths from the runners, or record them in a sink that every writer reports to, and build `output_files` from that record. Delete the folder scan.

## Validation

- A run into a folder that holds a GTR file from an earlier run, without a fitted model, does not list the GTR file
- An ancestral run with `--cdses=S` into a folder that holds the FASTA of CDS `M` lists only the FASTA of `S`
