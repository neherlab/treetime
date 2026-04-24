# Gap handling not implemented

v1 parses the `--keep-overhangs` CLI flag in the `ancestral` and `timetree` commands but does not wire it to any runtime behavior. Terminal gaps in input sequences are never replaced with the ambiguous character. This means v1 always behaves as if `--keep-overhangs` is set, while v0 defaults to filling terminal gaps (`fill_overhangs=True`).

The default v0 behavior (filling terminal gaps) is the scientifically correct default for most datasets: terminal gaps in aligned sequences typically represent missing data (incomplete sequencing), not true deletions. Treating them as informative gap characters biases ancestral reconstruction, GTR inference, and branch length estimation.

## Background

### What v0 does

v0's `seq2array()` ([packages/legacy/treetime/treetime/seq_utils.py#L152](../../packages/legacy/treetime/treetime/seq_utils.py#L152)) accepts a `fill_overhangs` parameter. When `True` (default), it finds the first and last non-gap positions in each sequence and replaces all gap characters outside that range with the ambiguous character (`N` for nucleotides, `X` for amino acids):

```python
if fill_overhangs:
    gaps = np.where(seq_array != '-')[0]
    if len(gaps):
        seq_array[:gaps[0]] = ambiguous      # leading gaps -> N
        seq_array[gaps[-1] + 1:] = ambiguous  # trailing gaps -> N
    else:
        seq_array[:] = ambiguous              # all-gap -> all N
```

This transformation happens at the I/O boundary, before any tree processing. The v0 CLI exposes `--keep-overhangs` which sets `fill_overhangs=False`. The wrapper code binds this as `fill_overhangs=not params.keep_overhangs` ([packages/legacy/treetime/treetime/wrappers.py#L333](../../packages/legacy/treetime/treetime/wrappers.py#L333)).

### What v1 does

v1 reads sequences verbatim from FASTA ([packages/treetime-io/src/fasta.rs#L101](../../packages/treetime-io/src/fasta.rs#L101)) without any terminal gap transformation. The `keep_overhangs` field is declared in:

- [packages/treetime/src/commands/ancestral/args.rs#L89](../../packages/treetime/src/commands/ancestral/args.rs#L89)
- [packages/treetime/src/commands/timetree/args.rs#L282](../../packages/treetime/src/commands/timetree/args.rs#L282)

Neither field is read by any production code. The destructuring in `run_ancestral_reconstruction()` ([packages/treetime/src/commands/ancestral/run.rs#L37](../../packages/treetime/src/commands/ancestral/run.rs#L37)) omits `keep_overhangs` via `..`.

### How gaps flow through v1 today

Gap characters (`-`) in input sequences are detected and tracked at multiple levels:

1. **Leaf payload construction.** `SparseNodePartition::new()` ([packages/treetime/src/representation/payload/sparse.rs#L43](../../packages/treetime/src/representation/payload/sparse.rs#L43)) calls `find_letter_ranges()` ([packages/treetime/src/seq/find_char_ranges.rs#L34](../../packages/treetime/src/seq/find_char_ranges.rs#L34)) to build `gaps: Vec<(usize, usize)>` and a corresponding `unknown` range list. The union forms `non_char` - positions excluded from the substitution model. `DenseNodePartition::new()` ([packages/treetime/src/representation/payload/dense.rs#L22](../../packages/treetime/src/representation/payload/dense.rs#L22)) similarly builds gap ranges.

2. **Fitch backward pass.** `run_fitch_backward()` ([packages/treetime/src/commands/ancestral/fitch.rs#L97](../../packages/treetime/src/commands/ancestral/fitch.rs#L97)) intersects child gap ranges to determine internal node gaps. Indel detection ([packages/treetime/src/commands/ancestral/fitch.rs#L212](../../packages/treetime/src/commands/ancestral/fitch.rs#L212)) compares child gap ranges against parent gap ranges to identify insertions and deletions.

3. **Marginal inference.** Both dense ([packages/treetime/src/representation/partition/marginal_dense.rs](../../packages/treetime/src/representation/partition/marginal_dense.rs)) and sparse ([packages/treetime/src/representation/partition/marginal_sparse.rs](../../packages/treetime/src/representation/partition/marginal_sparse.rs)) partitions use gap/non_char ranges to compute `edge_effective_length()` and to skip gap positions in `edge_subs()`.

4. **Profile mapping.** The alphabet maps gap to a uniform profile (identical to unknown) in `AlphabetConfig::create_profile_map()` ([packages/treetime/src/alphabet/alphabet_config.rs#L67](../../packages/treetime/src/alphabet/alphabet_config.rs#L67)), so gap positions carry no nucleotide signal during marginal inference. But gap ranges and unknown ranges are tracked separately and have distinct effects on indel detection and composition counting.

5. **Sequence composition.** `Composition` ([packages/treetime/src/seq/composition.rs](../../packages/treetime/src/seq/composition.rs)) counts character frequencies excluding the gap character. Terminal gaps that should be unknowns are currently counted as gaps, distorting the gap/unknown balance and the nucleotide frequency estimates used for GTR inference.

6. **Output.** Reconstructed ancestral sequences ([packages/treetime/src/representation/partition/marginal_dense.rs#L361](../../packages/treetime/src/representation/partition/marginal_dense.rs#L361)) overlay gap ranges back onto the argmax sequence. Terminal gaps that should be unknowns appear as `-` in output FASTA.

## Impact

The absence of terminal gap filling affects correctness and v0 parity:

- **Indel detection.** Terminal gaps in leaves are treated as true deletions during the Fitch backward pass. If a leaf has a 50-nucleotide terminal gap (incomplete sequencing), it generates a spurious large deletion indel on the edge to its parent. In v0 with default settings, these positions are `N` and do not trigger indel detection.

- **Effective sequence length.** `edge_effective_length()` subtracts gap/non_char positions from the total length. Terminal gaps in one leaf reduce the effective length for that leaf's edge, lowering the denominator for branch length estimation. With filling, these positions become unknowns, which have the same effect on effective length but a different effect on composition.

- **GTR inference.** Nucleotide frequency estimation via `Composition` excludes gap characters but counts canonical characters. Sequences with long terminal gaps have artificially short "informative" regions, skewing the base composition toward the non-terminal portion. Filling terminal gaps with `N` would exclude the entire terminal region uniformly.

- **Golden master parity.** All v0 comparison tests run with v0's default (`fill_overhangs=True`). v1 currently operates without filling, meaning numerical results diverge on datasets with terminal gaps.

## Proposed solution

Replace the boolean `keep_overhangs` flag with a three-mode enum and implement the terminal gap filling transformation at the I/O boundary.

### Gap handling modes

| Mode          | CLI value       | Behavior                                                                 | v0 equivalent                               |
| ------------- | --------------- | ------------------------------------------------------------------------ | ------------------------------------------- |
| Fill terminal | `only-terminal` | Replace leading and trailing gap characters with the ambiguous character | `fill_overhangs=True` (default)             |
| Fill all      | `all`           | Replace all gap characters with the ambiguous character                  | No v0 equivalent                            |
| None          | `none`          | Leave all gap characters unchanged                                       | `fill_overhangs=False` (`--keep-overhangs`) |

Default: `only-terminal` (matching v0 default).

`all` is a new mode not present in v0. It is useful for datasets where all gaps represent missing data rather than true deletions.

### Implementation approach

The transformation operates on `FastaRecord.seq` ([packages/treetime-io/src/fasta.rs#L16](../../packages/treetime-io/src/fasta.rs#L16)) after FASTA loading and before tree attachment. All downstream code remains unchanged - it sees sequences with gaps already replaced by unknowns where appropriate.

**Direct changes required:**

1. **Enum definition.** `GapFill { OnlyTerminal, All, None }` with `#[default] OnlyTerminal`, `Copy`, `Clone`, `ValueEnum`. Location: either `alphabet/` module or a shared args module.

2. **Transformation function.** `fn apply_gap_fill(seq: &mut Seq, mode: GapFill, gap: AsciiChar, unknown: AsciiChar)`. For `OnlyTerminal`: find first and last non-gap positions, fill the rest with `unknown`. For `All`: replace all `gap` with `unknown`. For `None`: no-op.

3. **CLI args.** Add `gap_fill: GapFill` to `ancestral/args.rs` and `timetree/args.rs`. Keep `--keep-overhangs` as a hidden flag (`#[clap(long, hide = true)]`) for v0 compatibility. Resolve early: if `--keep-overhangs` is set and `--gap-fill` is not explicitly provided, convert to `GapFill::None`. If both are provided, error.

4. **Wiring.** Call `apply_gap_fill()` on each `FastaRecord.seq` after loading in `ancestral/run.rs` (after [line 69](../../packages/treetime/src/commands/ancestral/run.rs#L69)) and `timetree/initialization.rs` (after FASTA loading).

**No changes required in downstream code.** The following components consume `Seq` data after loading and adapt automatically to changed gap/unknown distributions:

- `SparseNodePartition::new()`, `DenseNodePartition::new()` - gap range detection
- Fitch backward/forward passes - gap intersection, indel detection
- Marginal dense/sparse passes - gap-based position skipping
- `edge_effective_length()`, `edge_subs()` - gap/non_char union
- `Composition` - character counting
- GTR inference - base frequencies
- Sequence output - gap overlay

## Open questions

1. **Enum location.** `keep_overhangs` is duplicated across `ancestral/args.rs` and `timetree/args.rs`. Recommendation: shared `GapFill` enum in the `alphabet` module (it describes a property of the alphabet/sequence relationship), imported by command args.

2. **Command coverage.** `ancestral` and `timetree` load FASTA sequences. `optimize` also loads sequences. `clock` and `mugration` do not. Check whether `optimize` needs `--gap-fill` and add the flag if it loads FASTA.

### Testing

- Unit tests for `apply_gap_fill()`: each mode against sequences with leading gaps, trailing gaps, both, internal gaps, all-gap, no-gap, single character
- v0 parity: `apply_gap_fill(OnlyTerminal)` output must match v0 `seq2array(fill_overhangs=True)` for the same input
- Hidden flag resolution: `--keep-overhangs` alone resolves to `GapFill::None`, both flags together produces an error
- Golden master comparison suite before and after to quantify parity change on datasets with terminal gaps
