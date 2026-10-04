# `optimize` without `--alignment` fails with a missing-leaf-sequence error

`optimize` needs sequences, but `--alignment` is optional for every command. Without it, `optimize` reads zero records [packages/app-commands/src/commands/optimize/run.rs#L41](../../packages/app-commands/src/commands/optimize/run.rs#L41) and fails later in the Fitch pass:

```
Leaf sequence not found after alignment completion: 'EM_079497'
```

The message names the first leaf and does not say that the alignment is missing. `ancestral` rejects an empty `--alignment` list before reading with `--alignment is required: pass one or more FASTA files, or '-' to read standard input`.

## Expected behavior

`optimize` rejects an empty `--alignment` list before it reads any input, with the same message as `ancestral`.

## Validation

- A command test runs `optimize` with a tree and no alignment and checks the error message.
