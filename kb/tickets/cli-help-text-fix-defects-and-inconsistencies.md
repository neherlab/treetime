# Fix CLI help text defects and inconsistencies

Fix incorrect, misleading, and inconsistent `--help` output across all commands to serve human scientists, AI agents, and workflow engines.

## Tasks

### Defect fixes (high priority)

1. Remove `joint` from `--method-anc` possible values or hide it with a clap `hide` attribute and produce a parse-time error pointing users to `marginal`/`parsimony`
2. Use a different `--tree` description for commands where it is required (`optimize`, `ancestral`) -- remove the "If none is provided..." fallback clause. `prune` is already correct
3. Add a doc comment to `--prune-short` in `clock/args.rs` explaining what it does, or remove it if it should not exist on `clock`
4. Fix `--keep_root` underscore typo in `clock` command description to `--keep-root`
5. Replace `treetime/nuc_models.py` and `treetime/aa_models.py` references in `--model-params` with v1 Rust paths or documentation links
7. Update `ancestral` command description to not hardcode output filenames
8. Mark `homoplasy` and `arg` commands as unimplemented in their help text (add `[not yet implemented]` to the description)

### Inconsistency fixes

9. Unify `--clock-filter` description: pick one spelling of "interquartile", remove the redundant prose "Default=3.0"
10. Unify `--dense` descriptions across commands to a single shared doc comment
11. Unify `--vcf-reference` descriptions across commands
12. Unify `--max-iter` description capitalization ("Maximum" everywhere)
13. Unify `--branch-length-mode` description (remove stray "Note that" prefix in `clock`)
14. Standardize capitalization: sentence case for all flag descriptions

### UX improvements

15. Convert `--dense` from `true`/`false` string value to a boolean flag pair or tri-state with auto-detection
16. Replace LaTeX math notation in `--opt-method` possible values with plain text (`t`, `sqrt(t)`, `ln(t)`)
17. Add `--outdir` as a visible alias for `--output-all` for v0 migration
18. Add doc comment to `--pc` in `mugration` with a proper `#[default]` so `[default: 1.0]` appears in metadata
19. Add doc comment to `--detailed` in `homoplasy` explaining valid values
20. Restrict `--output-selection` possible values to the subset each command can produce

### Style fixes

21. Capitalize format names consistently: "Newick", "Nexus", "PHYLIP", "FASTA"
22. Fix grammar in `mugration --weights`: "probabilities of that" -> "probability that"
23. Add brief explanation for `--alphabet aa-no-stop`
24. Remove v0 internal references from `--smooth-initial-pi` and `--filter-uninformative-root` descriptions

### Workflow/automation fixes

25. Document that at least one output flag (`-O`/`--output-all` or a per-file flag) is required, or make it a clap validation
26. Generate per-command `--output-selection` possible values instead of the full superset
27. Document `--alignment` repeat syntax for multiple files (repeated flag, space-separated, etc.)
28. Make alignment stdin fallback text command-specific (remove from commands where stdin is irrelevant)
29. Replace root help docs link with v1 documentation or a disclaimer that linked docs cover v0
30. Add one minimal usage example per scientific command
31. Reformat `mugration --metadata` and `--weights` examples as structured multi-line help
32. Group root commands by purpose (Analysis, Tooling, Documentation)

### Grammar fixes

33. "If there's multiple input files" -> "If there are multiple input files" in `alignment.rs`
34. "gaussian" -> "Gaussian" in `timetree/args.rs`
35. "Higher numbers results" -> "Higher numbers result" in `mugration/args.rs`
36. Shorten `ancestral` root summary to one sentence, move detail to long help

## Scope

All changes are doc comments and clap attributes in `packages/treetime/src/commands/` and `packages/app-cli/src/cli/`. No runtime behavior changes except making `--method-anc=joint` a parse-time error instead of runtime error, fixing `--jobs` default display, and optionally adding clap validation for required output.

## Related issues

Source: [kb/issues/M-cli-help-text-defects.md](../issues/M-cli-help-text-defects.md)

Related:

- [kb/issues/N-timetree-dead-cli-flags.md](../issues/N-timetree-dead-cli-flags.md)
- [kb/issues/M-clock-dead-cli-arguments.md](../issues/M-clock-dead-cli-arguments.md)
- [kb/issues/H-timetree-tree-inference-unimplemented.md](../issues/H-timetree-tree-inference-unimplemented.md)
- [kb/issues/H-homoplasy-command-unimplemented.md](../issues/H-homoplasy-command-unimplemented.md)
- [kb/issues/N-timetree-polytomy-flags-no-conflict.md](../issues/N-timetree-polytomy-flags-no-conflict.md)
- [kb/issues/M-timetree-method-anc-ignored.md](../issues/M-timetree-method-anc-ignored.md)
