# Refactor Fitch compression long functions

Two functions in [packages/treetime/src/commands/ancestral/fitch.rs](../../packages/treetime/src/commands/ancestral/fitch.rs) exceed 170 lines each and handle multiple distinct phases inline. Both contain comments marking extractable blocks ("could be a function", numbered steps).

**Partial progress:** The indel processing logic (3-step backward and 3-loop forward) has been extracted into shared functions in [packages/treetime/src/commands/ancestral/fitch_indel.rs](../../packages/treetime/src/commands/ancestral/fitch_indel.rs) (`resolve_indels_backward`, `resolve_indels_forward`). Both `run_fitch_backward` and `run_fitch_forward` now delegate to these shared functions. The substitution and composition blocks remain inline.

## Affected functions

### `fn run_fitch_backward()` (174 lines, [L97-L270](../../packages/treetime/src/commands/ancestral/fitch.rs#L97-L270))

Logical blocks:

- Compute shared ranges (gaps, unknown, non_char, non_gap) from children
- Initialize parent sequence with FILL_CHAR / NON_CHAR
- Process variable positions (Fitch parsimony: intersection/union)
- Process fixed positions (children fixed, parent undecided)
- Process indels (3-step: child gap vs parent gap, variable-in-child propagation, all-children-gap collapse)
- Assemble `SparseNodePartition`

### `fn run_fitch_forward()` (191 lines, [L292-L482](../../packages/treetime/src/commands/ancestral/fitch.rs#L292-L482))

Logical blocks:

- Root branch: resolve variable states and indels by majority rule
- Non-root: copy parent to non_char positions, resolve variable positions (pick state or mutation), detect parent-only variable mutations
- Non-root: resolve indels (3 loops: variable_indel with parent tiebreak, consensus deletions, parent-only gaps)
- Non-root: adjust unknown composition
- All nodes: fill gaps/unknown in sequence, finalize root composition

## Approaches

Two approaches, from shallow to deep. Both preserve behavior.

### Approach A: extract helper functions

Extract each logical block listed above into a named function. The partition-iteration loop and the interleaving of substitution/indel/composition mutation remain as-is. Each helper still receives mutable references to the same shared state (`sequence`, `composition`, `gaps`, etc.).

Reduces line count per function. Does not change the coupling between concerns: partition locking, algorithm logic, and data structure mutation stay entangled across the extracted helpers.

### Approach B: separate computation from mutation

Both functions interleave three concerns: partition locking/iteration, Fitch algorithm logic (parsimony resolution), and data structure mutation (building `SparseNodePartition`, pushing subs/indels, updating composition). This approach separates them.

#### Split substitution resolution from indel resolution

Substitutions and indels are independent concerns that share mutable state only because they execute in the same scope. Substitutions operate on individual positions (state sets), indels operate on gap ranges. They don't depend on each other within a single pass step. Both update `composition` and `sequence`, which could be merged after both complete.

#### Pure computation -> apply pattern

Instead of `remove(&node.key)` + inline mutation + `insert(node.key, ...)`, compute a result struct from immutable inputs, then write once:

- Backward: `BackwardResult { gaps, unknown, non_char, sequence, variable, variable_indel }`
- Forward: `ForwardResult { resolved_sequence, subs, indels, composition }`

This makes the algorithm unit-testable without partition infrastructure and eliminates long mutable-borrow-through-locked-partition spans.

#### Tradeoff vs Approach A

Requires defining intermediate result types and splitting the composition updates (currently done incrementally during resolution) into a final merge step. More types, but each piece becomes independently testable and lock hold time shrinks.
