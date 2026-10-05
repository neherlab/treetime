# Homoplasy: mutation mapping and counting rules

`treetime homoplasy` reports how often each mutation and each site recurs across the branches of a tree. Recurrent mutations (homoplasies) point to selection, recombination, contamination, or fast-evolving sites. v1 keeps the statistics of v0 `scan_homoplasies()` ([`packages/legacy/treetime/treetime/wrappers.py#L82-L317`](../../packages/legacy/treetime/treetime/wrappers.py#L82-L317)) and changes how mutations reach them and how a few results are presented. The user approved each divergence below.

v1 code: `fn run_homoplasy()` in [`packages/app-commands/src/commands/homoplasy/run.rs`](../../packages/app-commands/src/commands/homoplasy/run.rs) maps the mutations; `fn run()` in [`packages/treetime/src/homoplasy/pipeline.rs`](../../packages/treetime/src/homoplasy/pipeline.rs) computes the statistics; [kb/algo/homoplasy.md](../algo/homoplasy.md) gives the formulas.

## Mutation mapping by `--method-anc`

v0 places mutations with joint maximum likelihood reconstruction (`infer_ancestral_sequences('ml', marginal=False)`). v1 has no joint reconstruction ([ancestral-joint-reconstruction-removed.md](ancestral-joint-reconstruction-removed.md)). `homoplasy` takes the `--method-anc` flag of `ancestral` (`marginal` by default, or `parsimony`) and runs the same reconstruction pipeline, so the branch mutations equal those of `ancestral` with the same flags.

The statistics need only the list of substitutions per branch, not the full joint configuration. On `data/zika/86`, v0 marginal reconstruction places the same 674 mutations on the same branches as v0 joint reconstruction, and v0 Fitch moves 8 of them. The golden-master test `test_gm_homoplasy_zika_86_matches_v0` compares v1 marginal with v0 joint on this dataset: every count, histogram, expected count, the log-likelihood difference, and every ranked list agree.

## Three mutation classes

v0 counts every parent-to-child difference except those involving `-` or `N`, so a change between a determined base and a partial ambiguity code counts as a substitution. On `data/flu/h3n2/20`, v0 lists `G451R` on two terminal branches as a homoplasy, although `R` (A or G) agrees with `G`.

v1 sorts every branch mutation into one of three classes with `fn classify_mutation()` ([`packages/treetime/src/homoplasy/classify.rs`](../../packages/treetime/src/homoplasy/classify.rs)):

- Substitution: both characters are canonical states (`A`, `C`, `G`, `T` for nucleotides). Only this class enters the v0 statistics: multiplicities, site hits, Poisson comparison, ranked lists, and homoplasic mutations per taxon
- Ambiguous change: a substitution with at least one character outside the canonical states, `N` and codes such as `R` or `Y` included. This class gets its own multiplicity histogram, ranked list, list of sites, and count per taxon
- Insertion or deletion: identified by kind, alignment columns, and the inserted or deleted characters. This class gets its own multiplicity histograms and ranked lists for all and for terminal branches, and the count per taxon of terminal indels that occur on two or more branches

The substitution class reads the branch mutations after unknown-state bridging (`struct UnknownMutationFilter` in [`packages/app-output/src/mutation_filter.rs`](../../packages/app-output/src/mutation_filter.rs)): `A` to `N` on one branch and `N` to `G` further down count once, as `A` to `G` on the lower branch. The ambiguous class reads the branch mutations before bridging, with leaf changes to `N` kept.

Gap filling turns the leading and trailing gaps of each sequence into `N` before reconstruction (`--gap-fill=only-terminal`, the v0 overhang filling), so these filled columns count as ambiguous changes on terminal branches. [kb/issues/N-homoplasy-filled-overhangs-dominate-ambiguous-changes.md](../issues/N-homoplasy-filled-overhangs-dominate-ambiguous-changes.md) tracks the effect on the ambiguous lists.

The "Sites mutated on several branches" tile of the `ancestral` results in the web and desktop apps (`fn ancestral_results()` in [`packages/app-commands/src/results/mutations.rs`](../../packages/app-commands/src/results/mutations.rs)) uses the same rule: it counts only substitutions between canonical states, so `G451R` and indel strings do not make a site recurrent there either.

## Smaller divergences

- Zero-hit count: v1 counts the sites without substitutions as the number of sites minus the sites with at least one substitution, in both indexing modes. v0 gets this count wrong when positions count from 1 ([kb/v0-errata/homoplasy-zero-hit-count-wrong-in-one-based-mode.md](../v0-errata/homoplasy-zero-hit-count-wrong-in-one-based-mode.md))
- DRM annotation: a substitution at a position of the `--drms` table gets the gene and the drug, and the amino-acid substitution only when the table lists its derived base. v0 crashes on an unlisted base ([kb/v0-errata/homoplasy-drm-annotation-crashes-on-unlisted-base.md](../v0-errata/homoplasy-drm-annotation-crashes-on-unlisted-base.md))
- List headers name the number of rows that `-n` selects ("The 20 most homoplasic mutations are:"). v0 always prints "ten", and writes "mutation" in the header of the terminal list
- Total tree length: v1 sums the branch lengths of the tree. v0 sets the root branch length to 0.001 and includes it, so the v1 total is 0.001 lower
- Poisson terms: v1 evaluates $\ln p_k$ with the log probability mass function, so a term with an underflowing $p_k$ stays finite and $p_k \ln p_k$ is 0 when $p_k$ is 0. v0 adds $10^{-100}$ to $p_k$ before the logarithm. The results agree unless $p_k$ underflows
- Tie order: mutations with the same multiplicity sort by position, then by ancestral and derived character. v0 keeps the insertion order of its dictionary at one position, which follows the tree traversal. Taxa with the same number of homoplasic mutations sort by name; v0 keeps the order of its tree
- Per-taxon table: with `--detailed`, the table has two more columns (ambiguous changes and recurrent indels on the terminal branch) and lists taxa that have one of the three counts above zero. When at least `-n` taxa carry homoplasic substitutions, the listed rows are the same as in v0

## Outputs

v0 prints the report to standard output and writes no files. v1 writes through the shared output plan of tree-writing commands: the statistics JSON (`homoplasy-stats`, `homoplasy.stats.json`), the text report (`homoplasy-report`, `homoplasy.report.txt`), and every tree format with the bridged branch mutations. `--output-all` writes the Newick and Nexus trees and both homoplasy files by default. The run also logs the report at info level (`-v`), so standard output stays free for outputs written to `-`. The statistics JSON contains the complete ranked lists, the terminal and per-taxon data, and the DRM annotations, whatever `-n` and `--detailed` select; `struct HomoplasyResult` in [`packages/app-commands/src/commands/homoplasy/result.rs`](../../packages/app-commands/src/commands/homoplasy/result.rs) documents its schema.

v0 flags that v1 drops: `--rng-seed` (replaced by `--seed`), `--verbose` (global `-v`), `--outdir` (v0 writes no files), and `--vcf-reference` (VCF input is unimplemented, [kb/issues/M-io-vcf-input-output-unimplemented.md](../issues/M-io-vcf-input-output-unimplemented.md)). `--tree` is required: v0 builds a tree with an external program when it is missing, and v1 has no tree building ([kb/issues/H-timetree-tree-inference-unimplemented.md](../issues/H-timetree-tree-inference-unimplemented.md)).
