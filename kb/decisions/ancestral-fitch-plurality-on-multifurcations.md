# Fitch reconstruction uses the plurality recurrence on multifurcations

v1 resolves a node's parsimony state set by intersecting its children's state sets and, when that
intersection is empty, retaining the states held by the greatest number of children. v0 retains
their union instead. On nodes with three or more children the union can keep states that are not of
minimum cost, so v1 produces different — and strictly more parsimonious — ancestral sequences and
mutation placements.

**Type**: Intentional scientific-correctness divergence from v0. Approved.

## Why v0's recurrence is wrong here

Let `g_v(x)` be the minimum number of changes in the subtree at `v` given state `x` at `v`, let
`M_i = min_y g_{c_i}(y)` be child `i`'s own minimum, and let `S_i` be its optimal set. A child
contributes `M_i` when `x ∈ S_i` and exactly `M_i + 1` otherwise, so

```text
g_v(x) = Σ_i M_i + k − count(x),   count(x) = #{ i : x ∈ S_i }
```

Minimising `g_v` therefore maximises `count`, and the node's minimum-cost set is
`{ x : count(x) = max }`.

Fitch's intersect-or-unite recurrence coincides with this for one or two children — a non-empty
intersection has `count = k`, and an empty one forces `max count = 1`, making the union the
plurality — but not beyond. For children `{C}`, `{C}`, `{A}` the union is `{A,C}`, which admits `A`
at a cost of two changes where `C` costs one.

The defect surfaces at a root polytomy because no parent state constrains the choice there. At
non-root nodes the forward pass prefers the parent's state whenever it is in the node's set, which
masks the error whenever the parent happens to be right.

## Measured effect

A 200-leaf root polytomy over real Ebola sequences (197,209 sites), scored over the 43,885
positions where every leaf is canonical:

| | substitutions |
| --- | --- |
| v0 recurrence (union + lowest byte) | 15,477 |
| plurality | 1,065 |

1,065 is the true optimum, computed independently from the alignment (`Σ_pos k − max_count`), and
the committed root state is of minimum cost at every one of those positions. The two disagree at
only 77 positions, but each costs up to `k − 1` extra changes, hence the 93% reduction. Across all
positions including gaps and ambiguity codes the totals are 58,930 against 4,364.

On a tree with a bifurcating root the output is unchanged: the 967-leaf dataset in `dbg/` produces
byte-identical mutation sets before and after, because its root has degree 2 and its internal
polytomies are masked by the parent-state preference.

## Performance

No measurable change. Alternating A/B over seven paired runs of `ancestral --method-anc parsimony`
on 967 leaves × 197,209 sites: 2.043 s against 2.046 s median, `+0.1%`.

The cost is dominated by the discovery scan, which is unchanged: comparing every child against
every position is `O(children · length)` and irreducible while nodes store dense sequences.
Instrumented, `fitch_backward` spent 374 ms of CPU in discovery against 8.7 ms in resolution, over
3,297 candidate positions across 471 internal nodes. Resolution is `O(children)` per candidate
rather than `O(1)` per disagreeing child, but candidates are far too few for that to matter.

## Implementation

The recurrence needs each child's optimal set plus a per-state count, not Sankoff cost vectors, so
the stored representation remains a single `StateSet` bitset and `SparseSeqInfo` is unchanged.
`StateSet::from_plurality` counts over the union of the input sets, whose width is bounded by the
alphabet, avoiding both a fixed-width counter array and its per-call zeroing.

The backward pass was reordered into discovery followed by resolution:

- `discover_fixed_disagreements_backward` keeps the child-major loop over whole sequences, writes
  the agreed state where children agree, and flags a position with `VARIABLE_CHAR` on the first
  disagreement. It builds no state sets and performs no map lookups.
- `resolve_variable_positions_backward` then runs over the flagged positions merged with the
  children's `fitch.variable` keys, recomputing each from every child.

Previously these ran in the opposite order and both folded child states into the same `variable`
map, because the second had no `VARIABLE_CHAR` guard. That double-counted children at any position
both observed, which is idempotent under a union and wrong under a count.

Note that leaves do carry `fitch.variable` entries: `SparseNodePartition::new` populates them for
every position holding a true IUPAC ambiguity code. Only `N` and `-` are excluded, via
`unknown` / `gaps` / `non_char`.

## Tie-breaking

Where several states remain equally parsimonious, the choice is free with respect to cost. It is
resolved deterministically, for reproducibility and stable golden masters, by
`Alphabet::first_canonical`, which orders by the configured alphabet rather than by byte value.
`StateSet::get_one()` orders by byte and is retained for callers that want that; it is no longer
used to commit a state to a sequence. For `nuc` the two coincide (`A`, `C`, `G`, `T`), so this part
changes nothing there. For `aa` they differ: canonical order lists `*` last while its byte (42)
sorts below `A`, so `get_one()` would commit a stop codon.

The forward pass's preference for the parent's state when it lies in the node's set is unchanged.
Given a minimal state set that preference is provably optimal rather than heuristic — matching the
parent avoids a change on the edge while any other member costs exactly one more — and it is what
pushes mutations toward the tips.

## Consequences

- `--model infer` derives GTR parameters from Fitch parsimony substitution counts, so the inferred
  model shifts on trees with multifurcations and marginal output shifts with it.
- Sparse marginal initialisation seeds `state` from `fitch.chosen_state` with a
  `profile.get_one()` fallback, so both of its inputs change.
- `test_fitch_polytomy` was rebaselined. At position 4 its `CDE` polytomy has child states
  `{T}`, `{T}`, `{C}`; the union committed `C` for three changes, the plurality commits `T` for
  two, which is minimal given that `B`, `C`/`D` and `E` hold three distinct states in separate
  subtrees.

## Code references

- `packages/treetime-primitives/src/bitset128.rs`: `BitSet128::from_plurality`
- `packages/treetime/src/alphabet/alphabet.rs`: `Alphabet::first_canonical`
- `packages/treetime/src/ancestral/fitch_sub.rs`: `discover_fixed_disagreements_backward`,
  `resolve_variable_positions_backward`, `choose_state`
- `packages/treetime/src/ancestral/fitch.rs`: discovery-then-resolution ordering
- `packages/legacy/treetime/treetime/treeanc.py#L638-L655`: v0's intersect-or-unite recurrence

## Related

- Source issue: [kb/issues/M-ancestral-fitch-polytomy-recurrence-not-minimum.md](../issues/M-ancestral-fitch-polytomy-recurrence-not-minimum.md)
- Ticket: [kb/tickets/ancestral-fitch-plurality-recurrence-on-multifurcations.md](../tickets/ancestral-fitch-plurality-recurrence-on-multifurcations.md)
- Hartigan, John A. 1973. "Minimum mutation fits to a given tree." _Biometrics_ 29:53-65. https://doi.org/10.2307/2529676
