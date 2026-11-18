# Codon substitution models are not implemented

TreeTime has no codon substitution model. Consistent nucleotide and amino-acid outputs do not need one ([M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md](M-ancestral-aa-outputs-disagree-with-translated-nucleotide-outputs.md)), but a codon model knows the genetic code: no history passes through a stop codon, and ω separates synonymous from nonsynonymous rates. Codon models fit RNA-virus genes better than nucleotide models; no study shows better ancestral states at viral divergence. Evidence: section "Evidence on codon substitution models" of [kb/reports/codon-aware-ancestral-reconstruction.md](../reports/codon-aware-ancestral-reconstruction.md#evidence-on-codon-substitution-models).

## Fix direction

[kb/proposals/codon-aware-ancestral-reconstruction.md](../proposals/codon-aware-ancestral-reconstruction.md) (axes 7 to 11) recommends an opt-in codon partition kind:

- **Model**: MG94xREV with one global ω, CF3x4 frequencies, and maximum-likelihood ω, GTR exchangeabilities and relative rate $\mu_c$, started from parsimony counts. No MG94 with F61 frequencies, which equals MG94 with F3x4
- **Units**: the generator normalized to one expected nucleotide substitution per nucleotide site, branch length $t = 3 \mu_c b$, and 3 sites per codon in `total_length`, `one_mutation`, the clock and the partition-rate normalization
- **Partition kind**: a first-class codon partition beside nucleotide partitions, declared by a `cds` key in `--partition` entries, with each coding column in exactly one model; codon model names stay out of the global nucleotide model enum
- **Tip observations**: ambiguous triplets as indicator rows over compatible sense codons; columns where a tip carries a stop codon, the annotated terminal stop codon, and codons that share nucleotides with another CDS stay in the nucleotide partition; partial gaps as missing codons
- **Representation**: sparse only. Before a 61-state run, change the shared per-state profiles of each message, the per-site heap allocations, and the L by n by n array per edge of transition counting (`fn get_branch_mutation_matrix()` [packages/treetime/src/gtr/infer_gtr/common.rs#L250-L279](../../packages/treetime/src/gtr/infer_gtr/common.rs#L250-L279)), and test the resolved-site threshold of 1e-4 at 61 states
- **Default**: codon decoding of the nucleotide reconstruction stays the default until the benchmark of [N-simulation-substitution-models-nucleotide-only.md](N-simulation-substitution-models-nucleotide-only.md) shows a gain

The section "Pitfalls found in an unmerged prototype" of [kb/reports/codon-aware-ancestral-reconstruction.md](../reports/codon-aware-ancestral-reconstruction.md#pitfalls-found-in-an-unmerged-prototype) lists the errors to avoid (branches 3 times too short, leaf rows of 1/k, erased stop codons, rejected overlaps).

> [!IMPORTANT]
> **Decision required.** Model family and frequencies: MG94xREV with CF3x4 only, GY94 with F61, or both. Tip observation model: columns with tip stop codons in the nucleotide partition, only the stop cell missing (IQ-TREE), or the column dropped (PAML, HyPhy). Whether a codon model ever becomes the default depends on the benchmark.

> [!IMPORTANT]
> **Investigation required.** Sparse 61-state memory (estimated 0.2 to 4 GB) and time (extrapolated to about 6 single-thread minutes on `mpox/clade-ii/2000`) are not measured, and ω identifiability on near-clonal data is not quantified.

## Validation

- MG94xREV with ω = 1 and a test code without stop codons against three codon-position nucleotide partitions with shared exchangeabilities: equal log-likelihoods, and codon posteriors equal to products of per-position posteriors. This also checks the factor 3 of the units
- Enumeration over all internal codon assignments on a tree with 2 or 3 internal nodes
- Matrix properties: rows sum to 0, stationarity, detailed balance, one expected nucleotide change per site per unit time; CF3x4 reproduces the observed per-position frequencies
- Dense against sparse codon runs on small datasets
