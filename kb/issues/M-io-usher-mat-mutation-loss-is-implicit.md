# UShER MAT output drops unrepresentable mutations implicitly

## Problem

The UShER MAT adapter silently changes the requested dataset after logging at most one warning per category. It drops all indels, all amino-acid substitutions, and every substitution containing a non-ACGT state [packages/treetime-io/src/tree_ir/usher.rs#L54-L87](../../packages/treetime-io/src/tree_ir/usher.rs#L54-L87). On input, it takes only the first value from the protobuf's repeated `mut_nuc` field and casts the signed position directly to `usize` [packages/treetime-io/src/tree_ir/usher.rs#L123-L133](../../packages/treetime-io/src/tree_ir/usher.rs#L123-L133).

UShER's protobuf mutation message has a single numeric position plus nucleotide state fields and no insertion/deletion event type [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/parsimony.proto#L4-L10)]. Dropping general indels is therefore a real format limitation. Choosing the first `mut_nuc` value is a separate reader policy: UShER combines every repeated value into an ambiguity-state bitset [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/src/mutation_annotated_tree.cpp#L565-L581)]. Negative positions are another distinct case: UShER treats them as masked mutations, while the TreeTime cast turns them into large ordinary coordinates [[src](https://github.com/yatisht/usher/blob/ac9c982d937c3bc43e2c16a0b73d30cf2b937118/src/matOptimize/mutation_annotated_tree_load_store.cpp#L321-L334)].

The current warnings do not let callers detect data loss programmatically, and a file is still produced successfully.

## Potential solutions

### A1. Track selection

- **Reject mixed nucleotide/amino-acid TreeIR:** fail whenever any non-nucleotide track is present.
- **Define MAT as a nucleotide-track projection:** ignore amino-acid tracks by contract, because the target format has no gene-track model. This must be documented as selection rather than an incidental warning.

### A2. Unrepresentable nucleotide events

- **Warn and drop:** preserve current behavior. Successful output does not imply preservation.
- **Fail by default:** reject indels, ambiguous substitutions, and unsupported repeated-state input with an aggregate diagnostic.
- **Explicit lossy mode:** fail by default and allow a caller-selected policy that reports the number and categories of dropped events.

### A3. Repeated `mut_nuc` input

- **Take the first state:** deterministic but discards schema information.
- **Map the state set to an IUPAC ambiguity code:** preserves supported nucleotide ambiguity when the set has a standard code.
- **Reject multi-state values:** avoids inventing a semantic interpretation but cannot read valid ambiguous MAT data.

### A4. Masked mutation input

- **Reject masked events:** prevents coordinate corruption but cannot preserve valid MAT masking information.
- **Preserve a masked-event variant:** extend TreeIR with an event that carries masking without pretending it is a nucleotide substitution.
- **Explicitly drop under lossy mode:** keep topology while reporting masked-event loss through the same aggregate policy as other unsupported data.

## Recommendation

Treat MAT output as an explicit nucleotide-track projection, fail by default on unrepresentable nucleotide events, and aggregate every loss category before returning the error. Add a lossy mode only after its CLI/API surface and reporting contract are approved. On input, map valid nucleotide state sets to IUPAC ambiguity codes, reject invalid or empty sets, and preserve masked events as a distinct TreeIR event if TreeIR remains the accepted boundary.

The checked numeric position conversion belongs to [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md); this issue concerns representational loss after coordinates are valid.

The writer also reconstructs `ref_nuc` from the root sequence. That is insufficient for recurrent mutations: after `A→G` on an ancestor, a descendant `G→T` must retain reference nucleotide `A`, not the immediate parent state `G`. Non-ASCII root bytes are currently filtered out, shifting every later coordinate. These are preservation failures governed by the same explicit MAT loss policy.

## Ticket readiness

No implementation ticket is ready. The default loss policy, whether a lossy mode should exist, and whether TreeIR must preserve masked events are externally visible decisions. The recommendation above must be approved before an implementation ticket is created.

## Related issues

- [M-core-mutation-representation-and-format-projection-inconsistent.md](M-core-mutation-representation-and-format-projection-inconsistent.md)
- [N-io-tree-ir-architecture-unapproved.md](N-io-tree-ir-architecture-unapproved.md)
