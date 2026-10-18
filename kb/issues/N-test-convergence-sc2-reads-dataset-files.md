# Optimize convergence tests read dataset files at run time

Both tests in `packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs` read their tree and alignment from `data/` at run time. The project testing rules forbid file input in unit tests, and the other optimize, timetree, prune and ancestral pipeline tests embed their datasets in the test binary through `enum TestDataset` [packages/treetime/src/test_utils/alignment.rs](../../packages/treetime/src/test_utils/alignment.rs).

## Current state

- `fn test_convergence_sc2_sparse_converges_on_sc2_2844()` joins `CARGO_MANIFEST_DIR` with `data/sc2/2844/tree.nwk` and `data/sc2/2844/aln.fasta.xz` and reads them with `fn nwk_read_file()` and `fn fasta_read_file()` [packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs#L26-L39](../../packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs#L26-L39)
- `fn test_convergence_sc2_flu_h3n2_20_converges()` reads `data/flu/h3n2/20/` the same way [packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs#L92-L105](../../packages/treetime/src/optimize/__tests__/test_convergence_sc2.rs#L92-L105)
- The SC2 files are 148 KB (tree) and 177 KiB compressed, 81.2 MiB decompressed (alignment). `enum TestDataset` already embeds `data/flu/h3n2/20/`

## Fix direction

- Read the flu case through `TestDataset::Flu20`
- Embed the SC2 dataset compressed, as `enum TestDataset` does, and decompress it in memory, or replace it with a smaller dataset that still exercises the convergence behavior the test checks. The choice depends on the test time and memory with the decompressed 81.2 MiB alignment

## Related

- `packages/treetime/src/io/__tests__/test_nwk.rs` lists the `tree.nwk` files under `data/` at run time
