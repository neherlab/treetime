# The run path check rejects the translation templates of augur workflows

`ancestral` stops with exit code 2 before it reads any input when `--translations` and `--output-translations` are templates in one folder whose file names differ by a suffix:

```
Error: --translations 'tr/%GENE.fasta' can name the same file as the output --output-translations 'tr/%GENE_withInternalNodes.fasta': the run would overwrite its own input
```

This is the naming of the augur workflows of Nextstrain: the rsv and seasonal-flu workflows pass `translations/%GENE.fasta` and `translations/%GENE_withInternalNodes.fasta`, so every `ancestral` step with translations of these workflows fails on v1, while v0 runs them.

## Cause

`fn check_run_paths()` [packages/app-commands/src/commands/shared/run_paths.rs#L30](../../packages/app-commands/src/commands/shared/run_paths.rs#L30) turns each template into a pattern in which every placeholder matches any file name, independently of the other template (`fn template_pattern()` [packages/app-commands/src/commands/shared/run_paths.rs#L194](../../packages/app-commands/src/commands/shared/run_paths.rs#L194), `fn patterns_intersect()` [packages/app-commands/src/commands/shared/run_paths.rs#L218](../../packages/app-commands/src/commands/shared/run_paths.rs#L218)). The input pattern `tr/*.fasta` and the output pattern `tr/*_withInternalNodes.fasta` intersect, because the input placeholder can take the value `X_withInternalNodes`. The two paths name the same file only when the run has two CDSs named `X` and `X_withInternalNodes`, and both templates take the CDS names of the run. The test case `both_templates` [packages/app-commands/src/commands/shared/__tests__/test_run_paths.rs#L28](../../packages/app-commands/src/commands/shared/__tests__/test_run_paths.rs#L28) asserts the pattern rule.

## Reproduction

```bash
printf '((A:0.1,B:0.2):0.05,C:0.3);\n' > tree.nwk
printf '>A\nATGAAATAA\n>B\nATGAAGTAA\n>C\nATGCAATAA\n' > aln.fasta
printf '##gff-version 3\nref\t.\tCDS\t1\t9\t.\t+\t0\tName=G;gene=G\n' > genome.gff3
mkdir -p tr && printf '>A\nMK*\n>B\nMK*\n>C\nMQ*\n' > tr/G.fasta
treetime ancestral --tree tree.nwk --alignment aln.fasta --annotation genome.gff3 \
  --translations 'tr/%GENE.fasta' --output-translations 'tr/%GENE_withInternalNodes.fasta' --output-node-data out.json
```

The run stops with the error above. With `--output-translations 'out/%GENE_withInternalNodes.fasta'` it finishes and writes `out/G_withInternalNodes.fasta`.

`test_run_paths_accept_augur_translation_templates_in_one_folder` in [packages/app-commands/src/commands/shared/__tests__/test_run_paths.rs](../../packages/app-commands/src/commands/shared/__tests__/test_run_paths.rs) asserts that the check accepts the two templates; it fails today and is ignored until the fix.

## Fix direction

Compare the concrete paths of the CDSs of the run: substitute each CDS name of `--cds-selection`, or of the annotation without it, into both templates, and reject the run only when a written path equals a read path. The check then needs the CDS names, so it runs once the annotation is read. The case `both_templates` changes with it: two templates collide only for CDS names that make their paths equal.

## Validation

- The ignored test passes and drops its `#[ignore]`
- The reproduction finishes with both templates in `tr/`
- A run whose CDSs are `X` and `X_withInternalNodes` still stops with the error
