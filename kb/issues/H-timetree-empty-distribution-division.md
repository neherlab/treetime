# Timetree aborts on dengue/500 and lassa/L/200 dividing by an empty distribution

`treetime timetree` with default options fails on `data/dengue/500` in round 1 and on `data/lassa/L/200` in round 2:

```
Error:
   0: When running round 1
   1: Timetree inference failed
   2: Cannot divide by empty distribution
```

The `rust` branch fails on both datasets with the same error, so the defect predates the current development branch.

## Symptom and reproduction

```bash
treetime timetree --tree=data/dengue/500/tree.nwk --dates=data/dengue/500/metadata.tsv \
  --aln=data/dengue/500/aln.fasta.xz --name-column=genbank_accession --output-all=<dir>
treetime timetree --tree=data/lassa/L/200/tree.nwk --dates=data/lassa/L/200/metadata.tsv \
  --aln=data/lassa/L/200/aln.fasta.xz --name-column=accession --output-all=<dir>
```

These are the `timetree/dengue/500/basic` and `timetree/lassa/L/200/basic` cases of `dev/smoke` (full tier).

## Impact and scope

- The command produces no output for these datasets.
- Both cases are declared expected failures in `dev/smoke.toml`, so a fix shows up as an unexpected pass until the declaration is removed.

## Mechanism

The division fails when one operand of a message division has empty support. An empty posterior arises when the message from the rest of the tree and a node's date constraint have disjoint support, the condition described in [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md). The source of the empty operand has not been traced.

## Related issues

- [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md)
- [H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md)
