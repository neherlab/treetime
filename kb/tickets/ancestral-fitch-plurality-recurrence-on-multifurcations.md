# Use the plurality (minimum-cost) Fitch recurrence on multifurcations

The sparse Fitch backward pass discards how many children carry each state. At a node with three or more children whose state sets do not share a common member, the surviving set is the union of the children's states, and the forward pass then resolves it with `StateSet::get_one()`, which returns the lowest ASCII byte. The result is a non-minimum reconstruction whose bias is visible at a root polytomy, where no parent state constrains the choice.

Reproducer: a five-leaf star with `l1..l4 = CGTT` and `l5 = AACA` yields parsimony root `AACA` (16 changes) against marginal root `CGTT` (4 changes, the minimum). With an internal polytomy (`(out,(l1..l5)inner)`, `out = AACA`) the node `inner` resolves to `AACA` at cost 4 where the minimum is 2.

Approved as an intentional divergence from v0, which uses the same intersection-or-union recurrence.

## Correctness argument

Let `g_v(x)` be the minimum number of changes in the subtree at `v` given state `x` at `v`, let `M_i = min_y g_{c_i}(y)`, and let `S_i = argmin_y g_{c_i}(y)` be child `i`'s optimal set. Then

```text
g_v(x) = Σ_i min( g_{c_i}(x), min_{y≠x} g_{c_i}(y) + 1 )
```

Under unit costs with integer `g`, the inner minimum is `M_i` when `x ∈ S_i` and exactly `M_i + 1` when `x ∉ S_i`. Therefore

```text
g_v(x) = C + k − count(x),    count(x) = #{ i : x ∈ S_i },   C = Σ_i M_i
```

so `argmin_x g_v(x) = argmax_x count(x)` and the node's optimal set is

```text
V(v) = { x : count(x) = max_y count(y) }
```

Three consequences that bound the work:

1. Only each child's **optimal set** and a per-state **count** are needed. No Sankoff cost vectors. The per-node stored representation stays a single `StateSet`; `SparseSeqInfo` does not change.
2. The existing intersection branch is **already exact**: members of a non-empty intersection have `count = k`, the maximum attainable, so `V(v) = intersection`.
3. `k = 2` is **already exact**: an empty intersection implies `max count = 1`, so `V(v) = union`.

New behaviour is therefore required only when `k ≥ 3` and the intersection is empty. `k` is the number of children carrying a character state at that position, not the node's degree.

## Prerequisite: invert the pass order, split discovery from resolution

The count is only well defined if each child contributes exactly once. Today it does not, because two functions both write `variable` and both fold in child states:

- `resolve_variable_positions_backward` (`fitch_sub.rs:18-79`) iterates positions drawn from the children's `fitch.variable` keys and folds in **every** child's state at those positions, using `StateSet::from_char(child.sequence[pos])` for children that are not variable there.
- `resolve_fixed_positions_backward` (`fitch_sub.rs:86-110`) then rescans **all children over the full sequence** and ORs each canonical child state into the same `variable` map.

The second function has no `VARIABLE_CHAR` guard — its only skips are `parent_state == child_state` and `parent_state == NON_CHAR` — so a position already promoted by the first function falls through and has its children folded in a second time. Under a bitwise union that is idempotent and harmless; under a count it is not.

The double-processing is confined to `VARIABLE_CHAR` positions. In the first function's `Unambiguous` case the second cannot re-widen: the intersection is a subset of every child's set, so a child fixed at `s` forces the resolved state to `s`, and the child that was variable carries `VARIABLE_CHAR` in its own sequence and fails `is_canonical`.

Both discovery sources must survive:

- Positions where some child carries a true IUPAC ambiguity code. `SparseNodePartition::new` (`partition/sparse.rs:43-49`) populates `fitch.variable` for every position where `Alphabet::is_ambiguous` holds, so **leaves do have `fitch.variable` entries**. `N` and `-` are excluded — they go to `unknown` / `gaps` / `non_char`.
- Positions where children hold differing canonical fixed states and no child is variable. These appear in no child's `fitch.variable` and are discovered only by comparing child sequences.

### Do not merge the two functions into one pass

`O(L·k)` byte comparisons are irreducible while each node stores a dense `sequence: Seq` of length `L`; discovering which positions differ requires touching every byte of every child. Only a sparse per-node diff representation would beat that, which is out of scope here.

The existing split exists so that the `O(L·k)` work iterates vectors and performs no position lookups. A single merged pass would have to consult `child.fitch.variable` per child per position — a `BTreeMap` lookup with poor locality inside the hot loop. Do not do this. Instead invert the order and separate discovery from resolution:

**Pass A — discovery** (`resolve_fixed_positions_backward`, runs **first**). Keep the child-major vector iteration exactly as it is: two streams per pass, no map operations. Detecting a disagreement does not require all children simultaneously, only that some child differs, so the cheap incremental loop suffices. On disagreement, write `VARIABLE_CHAR` and move on. **Do not construct a `StateSet` and do not touch the `variable` map.** This removes the `variable.entry(pos).or_insert_with(...)` `BTreeMap` operations that currently sit inside the `O(L·k)` loop.

**Collect.** One linear `O(L)` scan gathers `VARIABLE_CHAR` positions into a sorted `Vec<usize>`. No map, and sorted, so it merges cleanly against the children's already-sorted `fitch.variable` keys.

**Pass B — resolution** (`resolve_variable_positions_backward`, runs **second**) over discovered positions ∪ children's `fitch.variable` keys. It already re-collects every child's state set from scratch at each position it visits, so it is authoritative and counts each child exactly once. The plurality rule goes here and nowhere else.

Pass B visits more positions than today, but those are exactly the positions Pass A previously did map work for. Treat the net cost as unknown and measure it (see Validation); do not assume a speedup.

### Landmine to remove

`or_insert_with(|| stateset! {*parent_state})` (`fitch_sub.rs:104`) builds a state set from the sentinel when `*parent_state == VARIABLE_CHAR`, inserting `~` (ASCII 126) as a phantom state. It is unreachable today only because the first function always writes `VARIABLE_CHAR` and inserts into `variable` in the same step. It becomes reachable as soon as Pass A writes `VARIABLE_CHAR` without inserting, and a phantom state would then be counted. Removing this line is part of the fix.

## Implementation

### 1. Plurality primitive

Add a function that takes a slice of `StateSet` and returns the argmax-count set.

- Use a fixed-size counter array indexed by canonical alphabet index (`Alphabet::char`, `index_to_char`). Chosen over bit-sliced counting for readability and maintainability; revisit only if profiling shows it matters.
- No heap allocation on the hot path.
- Single pass to accumulate counts, single pass over the alphabet to select the maximum.
- Return the set of all states attaining the maximum, so downstream tie-break policy stays in one place.
- Assert or error on an empty input slice: at least one child is informative wherever the node is not itself `non_char`.

### 2. Pass A: discovery only

Reduce `resolve_fixed_positions_backward` to detection and reorder it to run before `resolve_variable_positions_backward`:

- Keep the child-major loop over `sequence.iter_mut().enumerate()` and the existing `parent_state == child_state` / `NON_CHAR` / `is_canonical(child_state)` conditions.
- `FILL_CHAR` still takes the child's state, as now.
- On any other disagreement, set `*parent_state = VARIABLE_CHAR` and nothing else. Drop the `variable.entry(...)` line entirely; the function no longer takes or returns a `variable` map.
- Add a `VARIABLE_CHAR` skip so an already-flagged position is not revisited by later children.

Then collect flagged positions with one linear scan of `sequence` for `VARIABLE_CHAR` into a sorted `Vec<usize>`.

### 3. Pass B: authoritative resolution

Extend `resolve_variable_positions_backward` to take the discovered positions in addition to the children's `fitch.variable` keys, and merge the two sorted sources into its position list. Per position it already does the right thing; only the empty-intersection branch changes:

1. Skip the position when the node is `non_char` there.
2. Collect each informative child's state set exactly once — `fitch.variable` entry when present, otherwise `StateSet::from_char(child.sequence[pos])` — skipping children that are `non_char` there, and preserving the existing `edge.transmission` skip.
3. Compute the intersection. If non-empty, use it (`Unambiguous` writes the state into `sequence`; `Ambiguous` records the set and writes `VARIABLE_CHAR`).
4. Empty intersection with `k == 2`: use the union (identical to the plurality result; keep as a fast path).
5. Empty intersection with `k ≥ 3`: use the plurality set.

Because Pass B recomputes from all children at every position it visits, it overwrites whatever Pass A left in `sequence`, so a flagged position gets exactly one count per child. Preserve the existing outputs unchanged: the returned `variable` map, `VARIABLE_CHAR` markers in `sequence`, and `FILL_CHAR` handling.

### 4. Tie-break policy

Apply one policy at all three sites that currently call `get_one()`:

- `resolve_root_forward` (`fitch_sub.rs:117-127`) — the root has no parent, so every member of `V(root)` is equally parsimonious.
- `resolve_nonroot_substitutions_forward` (`fitch_sub.rs:158` and `fitch_sub.rs:167`) — reached only when the parent state is not in the node's set.

Rules:

- **Keep the parent-state preference exactly as it is.** Given correct downpass sets, choosing the parent's state when it is in `V(v)` is provably optimal, not a heuristic: the child contributes `M_c` for `y = p ∈ S_c` and `M_c + 1` for any other choice. This is what pushes mutations toward the tips; do not alter it.
- **Remain deterministic.** Do not introduce randomised tie-breaking as the default; it would make output irreproducible between runs and destabilise golden masters without reducing cost.
- **Break remaining ties by canonical alphabet order, not ASCII order.** `StateSet::get_one()` is `bits.trailing_zeros()` (`bitset128.rs:155-168`) and never consults the alphabet. For `nuc` the two coincide (canonical `A,C,G,T`). For `aa` they do not: canonical lists `*` last but `*` is ASCII 42, below `A`, so any set containing a stop codon resolves to `*`. Either add an alphabet-aware selector for these call sites or give `Alphabet` a `first_canonical(StateSet)` helper. Leave `get_one()` itself alone if other callers depend on its current meaning.

## Validation

- Unit: the `C, C, A` three-child case returns `{C}`, not `{A, C}`.
- Unit: four children `{C}` and one `{A}` returns `{C}`; the reproducers above yield root `CGTT` and `inner` `CGTT`.
- Unit: genuine ties are preserved as sets — two children `{C}` and two `{A}` returns `{A, C}` — and the tie-break then selects deterministically.
- Unit: overlapping child sets are counted correctly — `k` children all `{A, C}` returns `{A, C}` with `count = k` for both, confirming that a simple majority over `k/2` does not imply uniqueness.
- Unit: a leaf carrying an IUPAC ambiguity code contributes its full state set once and not twice. Construct a position that both discovery sources would have found under the old code, so a double count would change the result.
- Unit: a position flagged by Pass A and also present in a child's `fitch.variable` appears exactly once in Pass B's merged position list.
- Unit: no `StateSet` ever contains a sentinel (`~`, `.`, or `' '`); assert over the returned `variable` map on inputs that promote positions in both passes.
- Benchmark: compare before and after on a tree with a large polytomy. The cost change is unknown by construction — Pass A loses its in-loop `BTreeMap` operations while Pass B gains positions — so record the measurement rather than assuming a direction. Measure by alternating the two binaries within one run; measuring the two sets minutes apart gives a spurious 30% regression from machine state alone.
- Unit: `k == 2` results are unchanged from the current implementation across intersecting and disjoint inputs.
- Unit: `aa` alphabet, a state set containing `*`, resolves to the lowest **canonical** state rather than `*`.
- Property: exhaustive brute-force minimum-change scoring over generated small multifurcations (degrees 3–6, alphabet size 4) confirms the returned set equals the true minimum-cost set at every node.
- Regression: confirm each new test fails against the current implementation before the change.

## Expected output changes

Recorded divergences to expect and approve, not defects:

- Golden-master ancestral outputs move on any tree containing a multifurcation.
- `--model infer` derives GTR parameters from Fitch parsimony substitution counts, so the inferred model changes and marginal output shifts with it. The probabilistic path will not stay byte-identical.
- Sparse marginal initialisation seeds `state` from `fitch.chosen_state` with a `profile.get_one()` fallback (`partition/marginal_passes.rs:283-293`); both inputs change.
- `fitch_cleanup` clears `variable` for internal nodes, so `run_fitch_reconstruction`'s `set_to_char` loop touches leaves only and is structurally unaffected.

Record the approved divergence in `kb/decisions/` once implemented, and un-ignore or re-baseline the affected golden masters.

## Status

Implemented. Outcome, measurements and the approved divergence are recorded in
[kb/decisions/ancestral-fitch-plurality-on-multifurcations.md](../decisions/ancestral-fitch-plurality-on-multifurcations.md).

## Related issues

- Source: [kb/issues/M-ancestral-fitch-polytomy-recurrence-not-minimum.md](../issues/M-ancestral-fitch-polytomy-recurrence-not-minimum.md) — note that this issue locates the defect only in `resolve_variable_positions_backward`; the polytomy-of-leaves path also runs through `resolve_fixed_positions_backward`, and both must change together.
