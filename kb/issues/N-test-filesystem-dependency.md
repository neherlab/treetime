# Tests with unnecessary filesystem dependency

Four unit tests read `data/` files when the underlying logic could use inline fixtures. This adds filesystem coupling, slows test execution, and makes tests fragile to working-directory changes.

All paths relative to `packages/treetime/src/`.

## `ancestral/__tests__/test_python_parity.rs` (1 test)

`test_root_sequence_matches_python_h3n2_na_20` reads `data/flu/h3n2/20/{tree.nwk,aln.fasta.xz}`. The tree is 1.5 KB (one line). The alignment is 20 sequences totaling 28 KB uncompressed FASTA, too large for a const literal but embeddable via `include_str!` with a decompressed fixture file.

## `clock/__tests__/test_clock_dengue100.rs` (2 tests)

Both tests read `data/dengue/100/{tree.nwk,metadata.tsv}`. The tree is 3.7 KB and the metadata is 101 lines (~2 KB). Embeddable as const, but the test's value is tied to this specific biological dataset (force_positive_rate behavior with all 198 root positions negative).

## `optimize/__tests__/test_convergence_sc2.rs` (1 test)

`test_convergence_sc2_flu_h3n2_20_converges` reads `data/flu/h3n2/20/{tree.nwk,aln.fasta.xz}`. Same dataset as python parity; shares the conversion approach.
