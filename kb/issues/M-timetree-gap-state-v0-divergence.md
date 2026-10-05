# Timetree treats the alignment gap as missing data, v0 as a fifth state

v0 timetree models the gap character `-` as a fifth nucleotide state for its default model and for the named models `JC69`, `F81` and `TN93`. v1 never gives the gap a state: the substitution model has the four nucleotide states, and gaps enter only as indel ranges. Sequence likelihoods, inferred GTR models and branch lengths in marginal branch-length mode can therefore differ from v0 on alignments with internal gaps. No decision records this difference.

## v0 behavior

- `fn create_gtr()` builds the model for `timetree` ([packages/legacy/treetime/treetime/wrappers.py#L375](../../packages/legacy/treetime/treetime/wrappers.py#L375)). The default `--gtr infer` ([packages/legacy/treetime/treetime/argument_parser.py#L194](../../packages/legacy/treetime/treetime/argument_parser.py#L194)) starts from `GTR.standard('jc', alphabet='nuc')` ([packages/legacy/treetime/treetime/wrappers.py#L55-L56](../../packages/legacy/treetime/treetime/wrappers.py#L55-L56))
- Named models take their alphabet from the model function defaults:
  - `JC69` and `F81` default to `alphabet='nuc'` ([packages/legacy/treetime/treetime/nuc_models.py#L18](../../packages/legacy/treetime/treetime/nuc_models.py#L18), [packages/legacy/treetime/treetime/nuc_models.py#L75](../../packages/legacy/treetime/treetime/nuc_models.py#L75))
  - `TN93` builds its GTR on `alphabets['nuc']` ([packages/legacy/treetime/treetime/nuc_models.py#L238](../../packages/legacy/treetime/treetime/nuc_models.py#L238)), although it sizes the equilibrium frequencies for four states ([packages/legacy/treetime/treetime/nuc_models.py#L230](../../packages/legacy/treetime/treetime/nuc_models.py#L230))
  - `K80`, `HKY85` and `T92` use `alphabets['nuc_nogap']` ([packages/legacy/treetime/treetime/nuc_models.py#L67-L70](../../packages/legacy/treetime/treetime/nuc_models.py#L67-L70), [packages/legacy/treetime/treetime/nuc_models.py#L146-L156](../../packages/legacy/treetime/treetime/nuc_models.py#L146-L156), [packages/legacy/treetime/treetime/nuc_models.py#L188](../../packages/legacy/treetime/treetime/nuc_models.py#L188))
- The `nuc` alphabet is `A, C, G, T, -` and maps `-` to a one-hot profile, so a gap is an observed state; `nuc_nogap` maps `-` to the all-ones profile, so a gap is missing data ([packages/legacy/treetime/treetime/seq_utils.py#L19-L54](../../packages/legacy/treetime/treetime/seq_utils.py#L19-L54))

## v1 behavior

- The nucleotide alphabet puts the gap in the undetermined set together with the unknown character ([packages/treetime/src/alphabet/alphabet.rs#L131](../../packages/treetime/src/alphabet/alphabet.rs#L131)), and `enum AlphabetName` has no variant with a gap state ([packages/treetime/src/alphabet/alphabet.rs#L304-L309](../../packages/treetime/src/alphabet/alphabet.rs#L304-L309))
- Gaps are tracked as indel ranges outside the substitution model; [kb/decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md](../decisions/ancestral-dense-gap-rule-observed-leaf-gaps.md) describes how internal nodes take their gap ranges from the leaves
- `--gap-fill` only replaces gaps by the unknown character (`fn apply_gap_fill()` in [packages/treetime/src/seq/gap_fill.rs](../../packages/treetime/src/seq/gap_fill.rs)), so no option restores the v0 model

## Scope

`ancestral`, `homoplasy`, and the unimplemented `arg` build their v0 model through the same `fn create_gtr()`, so the same difference applies to them. The decision for dense gap handling notes that v0 uses a different model, but no decision approves dropping the gap state.

## Open question

Record the four-state model with indel ranges as an intentional change in `kb/decisions/`, or add a five-state nucleotide alphabet for v0 parity.
