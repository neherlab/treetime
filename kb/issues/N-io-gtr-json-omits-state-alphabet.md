# The substitution model JSON of sequence models does not name its states

## Summary

The GTR output (`*.model.json`) is an object `{"models": [...]}` with one entry per partition that has a model, tagged by `kind`. A `trait` entry (mugration) carries `states`, the trait states in the order of `pi` and `W`. A `sequence` entry carries `partition`, `alignment`, `sites`, the optional `indel_rate` and `clock_rate`, and the model parameters `model_type`, `model_name`, `mu`, `pi`, `W` and `n_states`, but not the alphabet that orders `pi` and `W`. A reader cannot tell which frequency belongs to which nucleotide or amino acid without knowing the internal state order of the command that wrote the file.

## Evidence

- `struct SequenceModel` has no `states` field [`packages/treetime/src/gtr/get_gtr.rs#L43-L62`](../../packages/treetime/src/gtr/get_gtr.rs#L43-L62); `struct TraitModel` has one [`packages/treetime/src/gtr/get_gtr.rs#L64-L75`](../../packages/treetime/src/gtr/get_gtr.rs#L64-L75)
- `fn partition_models()` builds the sequence entries from `PartitionResult`, which holds the partition `alphabet` [`packages/app-commands/src/commands/shared/gtr_output.rs#L24-L40`](../../packages/app-commands/src/commands/shared/gtr_output.rs#L24-L40), [`packages/treetime/src/partition/result.rs#L17`](../../packages/treetime/src/partition/result.rs#L17)

## Impact

- The app shows only the model name and the overall rate `mu` of a sequence model, because it cannot label the frequencies (`struct SequenceModelSummary` [`packages/app-backend/src/results/models.rs#L40-L51`](../../packages/app-backend/src/results/models.rs#L40-L51)); trait models also show their states
- External tools that load the file must assume a state order
