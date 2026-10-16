# Sparse profiles allocate for alphabets above four states

Sparse marginal messages keep each site profile and each fixed profile inline in a `SmallVec` with the capacity `INLINE_STATES = 4` (`type SiteProfile` [packages/treetime/src/partition/storage/sparse.rs#L17-L19](../../packages/treetime/src/partition/storage/sparse.rs#L17-L19)), the canonical state count of the nucleotide alphabet. With an amino-acid alphabet (20 or 21 states), every profile exceeds the inline capacity and `SmallVec` moves it to the heap. Each variable site and each fixed state of a message then takes one allocation, and the list of fixed profiles (`struct FixedProfiles`, inline capacity of 4 entries) takes one more. Results are the same for every alphabet; only the allocation count differs.

A larger inline capacity is not a remedy: a capacity of 20 makes every nucleotide message 3 to 4 times larger.

> [!IMPORTANT]
> **Decision required.** Two designs keep the allocation count independent of the alphabet size:
>
> - One profile buffer per message, sized from the alphabet: the entries refer to rows of one `Vec<f64>`, so a message takes 2 to 3 allocations for any alphabet
> - Kernels generic over the state count, with the inline capacity chosen per alphabet where a command selects its alphabet

> [!IMPORTANT]
> **Investigation required.** The cost is not measured. The only protein alignments in the repository are the 20-sequence translations in `data/rsv/a/20/translations/`. Measure allocation counts and run time of sparse `ancestral` with the amino-acid alphabet on a larger protein alignment before choosing a design.

## Related issues

- [N-ancestral-marginal-array-kernels-allocate.md](N-ancestral-marginal-array-kernels-allocate.md): per-site temporaries of the marginal kernels
