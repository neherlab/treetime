# Input name matching and validation for tree-data reconciliation

Analysis commands reconcile names between independent input sources: tree leaf names (Newick), sequence names (FASTA), date labels (metadata TSV), and trait labels (metadata TSV). This reconciliation has correctness defects, poor diagnostics, and a performance regression relative to v0. This proposal introduces indexed lookup, batch error reporting, duplicate detection, and optional fuzzy suggestions as four independent improvements.

## Problem

### Defective call sites

Six call sites perform name reconciliation. Three have known defects:

| Call site                                                                                                    | Lookup                      | Batch errors           | Duplicate detection   |
| ------------------------------------------------------------------------------------------------------------ | --------------------------- | ---------------------- | --------------------- |
| `fitch.rs:54-101` [[src](../../packages/treetime/src/ancestral/fitch.rs#L54-L101)]                           | linear scan O(n\*m)         | no (first miss aborts) | no                    |
| `marginal_dense.rs:319-351` [[src](../../packages/treetime/src/partition/marginal_dense.rs#L319-L351)]       | linear scan O(n\*m)         | no (first miss aborts) | no                    |
| `attach.rs:28-89` [[src](../../packages/treetime/src/ancestral/attach.rs#L28-L89)]                           | `BTreeSet` O(n log n)       | warnings + count       | no                    |
| `date_constraints.rs:21-90` [[src](../../packages/treetime/src/clock/date_constraints.rs#L21-L90)]           | `BTreeMap::get` O(log n)    | warnings               | n/a (map keys unique) |
| `marginal_discrete.rs:214-247` [[src](../../packages/treetime/src/partition/marginal_discrete.rs#L214-L247)] | `IndexSet::difference` O(n) | bidirectional batch    | n/a (map keys unique) |
| `marginal_discrete.rs:47-90` [[src](../../packages/treetime/src/partition/marginal_discrete.rs#L47-L90)]     | `BTreeMap::get` O(log n)    | via validate above     | n/a                   |

The two sequence attachment functions (`fitch.rs` and `marginal_dense.rs`) use `aln.iter().find(|fasta| fasta.seq_name == leaf_name)` -- a linear scan per leaf node, yielding O(n\*m) total comparisons where n = leaves and m = sequences. They abort on the first unmatched name, forcing the user to fix one mismatch at a time and re-run. Neither detects duplicate names in the input, so duplicate tree leaves silently receive the first matching sequence and duplicate FASTA names silently shadow each other.

The mugration validator `validate_trait_names` (`marginal_discrete.rs:214-247` [[src](../../packages/treetime/src/partition/marginal_discrete.rs#L214-L247)]) already solves all three problems: it uses `IndexSet::difference` for O(n) bidirectional batch reporting. This is the internal model the sequence attachment sites should follow.

### Performance regression from v0

v0 stores the alignment as a Python `dict` keyed by name (`treeanc.py:395-444` [[src](../../packages/legacy/treetime/treetime/treeanc.py#L395-L444)]), giving O(1) lookup per leaf. v1's linear scan is a regression:

| Sequences | v1 current (linear scan) | v0 / indexed  |
| --------- | ------------------------ | ------------- |
| 1,000     | 1,000,000 comparisons    | 1,000 lookups |
| 10,000    | 100,000,000              | 10,000        |
| 100,000   | 10,000,000,000           | 100,000       |
| 1,000,000 | 1,000,000,000,000        | 1,000,000     |

Pandemic-scale datasets reach millions of sequences (SARS-CoV-2 phylogenetics with UShER processes 14M+ genomes). The quadratic loop becomes the dominant cost before any phylogenetic algorithm runs.

### Poor diagnostic quality

Current error on mismatch:

```
Leaf sequence not found after alignment completion: 'A/California/07/2009'
```

The user must manually diff thousands of names across tree and FASTA files to find whether the cause is a typo, case mismatch, underscore/space swap, truncation, or entirely wrong file. The error provides no list of near-matches, no count of total failures, and no indication of whether the problem is systematic or isolated.

## Ecosystem survey

Nine phylogenetics tools were inspected at the source level to understand how the ecosystem handles tree-alignment name reconciliation. For each tool, the actual comparison function, duplicate detection logic, and error reporting code were traced in source. Every claim below links to the specific source location.

### Matching semantics

All tools use exact string equality for the primary lookup. Where normalization exists, it is applied symmetrically to both sides as a preprocessing step before the exact match -- not as a fuzzy fallback. This distinction matters: normalization changes which data attaches to which tree node and therefore affects scientific results. Fuzzy suggestions (axis 3 in this proposal) are purely diagnostic and do not change matching behavior.

| Tool        | Comparison                                             | Normalization                             | Source                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| ----------- | ------------------------------------------------------ | ----------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| v0 TreeTime | `name not in dict` (exact)                             | none                                      | [[src](https://github.com/neherlab/treetime/blob/ea8a21a74d8cbfd36466caa89f60e4736621de82/packages/legacy/treetime/treetime/treeanc.py#L395-L444)]                                                                                                                                                                                                                                                                                                                                     |
| augur       | delegates to v0                                        | space-split for IQ-TREE compat only       | [[src](https://github.com/nextstrain/augur/blob/024292af6daf1f8aae7e8d690bd2f944e1c7322d/augur/align.py#L197-L211)]                                                                                                                                                                                                                                                                                                                                                                                                     |
| IQ-TREE 2   | exact after `renameString()`                           | non-alphanum -> `_` on both sides         | [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/utils/tools.cpp#L503-L511)], [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/tree/phylotree.cpp#L349-L383)]                                                                                                                                                          |
| BEAST2      | `String.equals()` (exact)                              | none                                      | [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/alignment/Alignment.java#L321)], [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/TreeParser.java#L424)]                           |
| BEAST v1    | `String.equals()` (exact)                              | FASTA header trim only                    | [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/util/TaxonList.java#L114)], [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evomodel/treelikelihood/BeagleTreeLikelihood.java#L490)] |
| MrBayes     | exact after `tolower()`                                | case fold on both sides (always)          | [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/utils.c#L1707-L1726)]                                                                                                                                                                                                                                                                                                                                                                                                       |
| DendroPy    | `str.lower()` when `is_case_sensitive=False` (default) | unquoted `_` -> space per Newick standard | [[src](https://github.com/jeetsukumaran/DendroPy/blob/75bc2aea3450ef48a7aba5017551318f6586b0f6/src/dendropy/datamodel/taxonmodel.py#L538)]                                                                                                                                                                                                                                                                                                                                        |
| gotree      | Go map lookup (exact)                                  | none                                      | [[src](https://github.com/evolbioinfo/gotree/blob/83cca67fdcf64af47f6f3318d3eebdcb417d45a4/tree/tree.go#L454-L468)]                                                                                                                                                                                                                                                                                                                                                                                                      |
| ETE 4       | Python dict lookup (exact)                             | none                                      | [[src](https://github.com/etetoolkit/ete/blob/9562dfb6a02795dfda3975be5025f83b8dc884b1/ete4/phylo/phylotree.py#L346-L369)]                                                                                                                                                                                                                                                                                                                                                                                   |

IQ-TREE's `renameString()` is the most aggressive normalizer: it replaces any character outside `[A-Za-z0-9_.-]` with `_`, applied to both tree and alignment names at parse time. MrBayes applies `tolower()` unconditionally through its central `StrCmpCaseInsensitive()` function used in all name lookups. DendroPy is the only tool implementing the Newick-standard underscore-to-space rule (`preserve_underscores=False` default in its Newick reader).

### Duplicate detection

| Tool        | Tree duplicates                   | Alignment duplicates                 | Source                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ----------- | --------------------------------- | ------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MrBayes     | hard error (all contexts)         | hard error                           | [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/command.c#L7234-L7285)]                                                                                                                                                                                                                                                                                                                                                                            |
| gotree      | hard error on index build         | n/a                                  | [[src](https://github.com/evolbioinfo/gotree/blob/83cca67fdcf64af47f6f3318d3eebdcb417d45a4/tree/tree.go#L454-L468)]                                                                                                                                                                                                                                                                                                                                                                                 |
| BEAST2      | hard error (first)                | hard error (first)                   | [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/TreeParser.java#L462-L471)], [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/alignment/Alignment.java#L321)] |
| IQ-TREE 2   | none                              | batch (sorted adjacent scan)         | [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/alignment/alignment.cpp#L136-L161)]                                                                                                                                                                                                                                                                                                                                                               |
| augur       | none                              | batch (sorted set, all dupes listed) | [[src](https://github.com/nextstrain/augur/blob/024292af6daf1f8aae7e8d690bd2f944e1c7322d/augur/io/sequences.py#L506-L540)]                                                                                                                                                                                                                                                                                                                                                              |
| BEAST v1    | silent dedup in `Taxa.addTaxon()` | none                                 | [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/util/TaxonList.java#L114)]                                                                                                                                                                                                                                                                                                                           |
| v0 TreeTime | none                              | silent (dict key overwrite)          | [[src](https://github.com/neherlab/treetime/blob/ea8a21a74d8cbfd36466caa89f60e4736621de82/packages/legacy/treetime/treetime/treeanc.py#L419)]                                                                                                                                                                                                                                                                                                                     |
| ETE 4       | error in EvolTree only            | silent auto-rename (`2_name`)        | [[src](https://github.com/etetoolkit/ete/blob/9562dfb6a02795dfda3975be5025f83b8dc884b1/ete4/parser/fasta.py#L8-L50)]                                                                                                                                                                                                                                                                                                                                                                             |

MrBayes and gotree are the strictest: any duplicate name is a hard error. BEAST2 detects duplicates in both sources but reports only the first. IQ-TREE and augur detect alignment duplicates in batch but not tree duplicates. v0, BEAST v1, and ETE handle duplicates silently (overwrite, dedup, or rename), masking data errors.

### Batch error reporting

| Tool        | Behavior                                                                                         | Source                                                                                                                                                                                                                                                                                                   |
| ----------- | ------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| IQ-TREE 2   | bidirectional batch: all alignment-not-in-tree + all tree-not-in-alignment reported before abort | [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/tree/phylotree.cpp#L349-L383)]                                                             |
| augur       | batch for FASTA duplicates (sorted list); one-at-a-time elsewhere                                | [[src](https://github.com/nextstrain/augur/blob/024292af6daf1f8aae7e8d690bd2f944e1c7322d/augur/io/sequences.py#L506-L540)]                                                  |
| BEAST2      | set-difference warnings (non-fatal, both directions); one-at-a-time for hard errors              | [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/TreeParser.java#L162-L181)] |
| v0 TreeTime | per-leaf warning + cumulative count + abort when >1/3 missing                                    | [[src](https://github.com/neherlab/treetime/blob/ea8a21a74d8cbfd36466caa89f60e4736621de82/packages/legacy/treetime/treetime/treeanc.py#L395-L444)]    |
| All others  | one-at-a-time (abort on first mismatch)                                                          | see matching semantics table                                                                                                                                                                                                                                                                             |

IQ-TREE's `setAlignment()` is the best implementation: it iterates both directions (alignment -> tree, tree -> alignment), collects all mismatches, and reports them in a single error. This is the external model for batch reporting.

### Fuzzy suggestions

No surveyed tool provides fuzzy matching suggestions for unmatched names. This would be novel in the phylogenetics ecosystem.

### Newick underscore-space standard

The Newick specification (Felsenstein, 1986) states that underscores in unquoted labels represent spaces: "it is assumed that an underscore character stands for a blank; any of these in a name will be converted to a blank when it is read in" ([source](https://phylipweb.github.io/phylip/newicktree.html)). See also [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md).

DendroPy and scikit-bio implement this conversion by default. All other surveyed tools (v0, v1, augur, BEAST, IQ-TREE, gotree, ETE, MrBayes) keep underscores verbatim. The v1 Newick parser does not convert underscores ([[src](https://github.com/neherlab/treetime/blob/ea8a21a74d8cbfd36466caa89f60e4736621de82/packages/util-newick/src/parse.rs)]).

## Design axes

Four independent decision axes. Each can be adopted, deferred, or rejected independently.

### Axis 1: Lookup data structure

Replace linear-scan `.find()` with indexed lookup. Pure performance fix, no behavior change.

The two affected sites are `fitch.rs:72-76` [[src](../../packages/treetime/src/ancestral/fitch.rs#L72-L76)] and `marginal_dense.rs:331-335` [[src](../../packages/treetime/src/partition/marginal_dense.rs#L331-L335)], both using `aln.iter().find(|fasta| fasta.seq_name == leaf_name)`.

| Sequences | Current (linear scan)         | HashMap           |
| --------- | ----------------------------- | ----------------- |
| 1,000     | 1,000,000 comparisons         | 1,000 lookups     |
| 100,000   | 10,000,000,000 comparisons    | 100,000 lookups   |
| 1,000,000 | 1,000,000,000,000 comparisons | 1,000,000 lookups |

#### Ordering impact of Vec -> HashMap

Replacing `Vec<FastaRecord>` with `HashMap<String, FastaRecord>` destroys insertion order. Code analysis confirms this is safe:

- All leaf-to-sequence matching is by name, not position (`.find()` by `seq_name`)
- `FastaRecord.index` (set during reading) is unused -- never read in production
- `complete_alignment_for_leaves` builds a `BTreeSet` from names internally, then appends synthetic records; append position is irrelevant since downstream lookup is by name
- `create_mask` ORs non-ambiguous positions across all records -- commutative, order-independent
- `get_common_length` checks all records have equal length -- set check, order-independent
- Output FASTA order is determined by tree traversal (preorder), not input FASTA order
- Node data JSON is keyed by `GraphNodeKey`, not by FASTA position

One test generator (`prop_generators/alignment.rs:94`) asserts positional ordering of generated records -- cosmetic, easy to adjust.

#### Where to build the index

The FASTA data flows through three stages:

1. `read_many_fasta()` (`treetime-io/src/fasta.rs` [[src](../../packages/treetime-io/src/fasta.rs)]) produces `Vec<FastaRecord>`
2. `complete_alignment_for_leaves()` (`attach.rs:28-89` [[src](../../packages/treetime/src/ancestral/attach.rs#L28-L89)]) checks for missing leaves (builds a temporary `BTreeSet` of names internally), appends synthetic ambiguous records for missing tips, returns `Vec<FastaRecord>`
3. `attach_seqs_to_graph()` / `attach_sequences()` iterates leaves, does the linear `.find()` lookup

##### Option A: Index at attachment boundary

Build `HashMap<&str, &FastaRecord>` inside `attach_seqs_to_graph` and `attach_sequences`, from the incoming `&[FastaRecord]`.

Pros:

- Smallest diff -- changes scoped to two functions
- No API change, all callers unaffected
- No risk to downstream consumers of `Vec<FastaRecord>`

Cons:

- `complete_alignment_for_leaves` retains its own redundant `BTreeSet` index
- FASTA duplicate detection happens late (at attachment, not at read time)
- Two attachment sites must each build the index independently (duplication)

##### Option B: Index at parse time

Have `read_many_fasta` return `HashMap<String, FastaRecord>`. The indexed representation propagates through the entire pipeline.

Pros:

- Single index, built once, reused everywhere
- `complete_alignment_for_leaves` uses the map directly -- no redundant `BTreeSet`
- FASTA duplicates detected at the earliest point (during reading)
- Cleaner data model: the map enforces name uniqueness structurally

Cons:

- `Vec<FastaRecord>` appears in ~50 files (mostly tests); type change propagates
- Tests constructing inline alignments need adjustment (build a `HashMap` instead of `vec![]`)

##### Recommendation: Option B

The type change propagation is mechanical (largely `vec![...]` -> map construction in tests). The benefits are structural: duplicate detection at the earliest point, one index instead of three redundant ones (`BTreeSet` in `complete_alignment_for_leaves`, `.find()` in `fitch.rs`, `.find()` in `marginal_dense.rs`), and the data model itself prevents the class of duplication bugs.

### Axis 2: Error collection and reporting

Replace one-at-a-time abort with batch collection and structured error messages. Diagnostic improvement only -- no change to matching semantics.

Pros:

- Users see all problems at once instead of fixing one, re-running, fixing the next
- Bidirectional reporting (IQ-TREE model) makes systematic causes visible (e.g. all names differ by underscore/space)
- Duplicate detection prevents silent wrong-sequence attachment -- a correctness defect, not just UX
- Both the external model (IQ-TREE) and internal model (mugration `validate_trait_names`) already demonstrate this approach

Cons:

- Error messages become longer for large mismatch counts; need truncation for readability
- Collecting all mismatches before aborting means running the full matching loop even when the first mismatch indicates total failure (wrong file). Mitigated by the >1/3 threshold from v0

Recommendation: implement both 2a and 2b. These are correctness fixes, not optional features.

#### 2a: Batch mismatch reporting

Collect all unmatched names, report in a single error. Bidirectional: both "leaves not in alignment" and "alignment names not in tree", following the patterns of IQ-TREE's `setAlignment()` ([[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/tree/phylotree.cpp#L349-L383)]) externally and mugration's `validate_trait_names` (`marginal_discrete.rs:214-247` [[src](../../packages/treetime/src/partition/marginal_discrete.rs#L214-L247)]) internally.

Example error output:

```
Name reconciliation failed: 3 tree leaves not in alignment, 2 alignment names not in tree.

Tree leaves not in alignment:
  'A/California/07/2009'
  'A/Texas/50/2012'
  'B/Brisbane/60/2008'

Alignment names not in tree:
  'A_California_07_2009'
  'A_Texas_50_2012'
```

When both sides have unmatched names showing the same transformation pattern (underscore/space, case), the pattern becomes visible in the output, letting the user diagnose systematic causes without manual diffing.

#### 2b: Duplicate detection

During index construction, detect and report duplicate names in both inputs before matching begins. Duplicates in either source are a data error: duplicate tree leaves cause ambiguous sequence attachment, and duplicate FASTA names cause silent shadowing.

Example error output:

```
Duplicate leaf names in tree: 'USA' (2 occurrences), 'CHN' (3 occurrences)
Duplicate sequence names in alignment: 'sample_1' (2 occurrences)
```

Detection is free with `HashMap`: on `insert`, check the return value for a displaced entry. Track collisions in a `Vec<(&str, usize)>` and report all before aborting.

### Axis 3: Fuzzy suggestions for unmatched names

For each unmatched name, compute string similarity against names from the other source and report top-N near-matches in the error message. Advisory only -- matching semantics remain exact. No surveyed tool provides this.

Pros:

- Turns opaque "name not found" errors into actionable diagnostics (shows what the name is close to)
- Zero risk to correctness -- suggestions are displayed in error output, never used for matching
- Zero new dependencies (`strsim` already in dependency tree via `clap`)
- Novel in the phylogenetics ecosystem -- no surveyed tool offers this

Cons:

- O(U \* M) search for U unmatched names against M candidates; requires capping and filtering for large datasets (see performance analysis below)
- Suggestions can be misleading when names are structurally similar but semantically distinct (`A/California/07/2009` vs `A/California/08/2009` are different strains, not typos)
- Adds complexity to error formatting code

Recommendation: implement with the capping strategy described below. The UX benefit is high and the risk is zero (advisory only).

#### Algorithm

Taxon names are short strings (10-100 characters, ASCII-dominated). Relevant similarity metrics:

| Metric              | Complexity per pair | Properties                                                                                                                                  |
| ------------------- | ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| Jaro-Winkler        | O(n)                | Prefix-weighted. Phylogenetic names share long prefixes (`A/California/07/2009` vs `A/California/7/2009`), making this metric a natural fit |
| Levenshtein         | O(n\*m)             | Edit operations (insert, delete, substitute). Intuitive interpretation                                                                      |
| Damerau-Levenshtein | O(n\*m)             | Adds transposition. Catches `AB` vs `BA` swaps                                                                                              |

Jaro-Winkler is the recommended primary metric: O(n) per pair, produces a 0-1 similarity score, and prefix weighting aligns with phylogenetic naming conventions where the informative variation is at the end of the name.

#### Dependency

`strsim 0.11.1` is already in the dependency tree as a transitive dependency of `clap` (via the `suggestions` feature, which `clap` uses for "did you mean" hints on mistyped subcommands). Using `strsim` directly adds zero new dependencies to the build. It provides Jaro-Winkler, Levenshtein, Damerau-Levenshtein, and normalized variants.

For higher-throughput needs if profiling identifies a bottleneck:

- `rapidfuzz` (Rust port): `BatchComparator` caches preprocessing for 1-to-many comparisons with `score_cutoff` for early termination
- `editdistancek`: bounded Levenshtein at O(d \* min(n,m)), sublinear when most pairs are distant
- `triple_accel`: SIMD-accelerated (AVX2/SSE4.1), `u8`-only, designed for bioinformatics sequence comparison

#### Performance at scale

Fuzzy search runs only for unmatched names, not for all names. Cost is O(U _ M _ c) where U = unmatched count, M = candidate count from the other source, c = per-pair cost.

| Scenario                  | U (unmatched) | M (candidates) | Pairs          | Per-pair cost | Estimated wall time |
| ------------------------- | ------------- | -------------- | -------------- | ------------- | ------------------- |
| Few typos, small dataset  | 5             | 1,000          | 5,000          | O(n), n~30    | <1ms                |
| Few typos, large dataset  | 5             | 100,000        | 500,000        | O(n), n~30    | ~5ms                |
| Few typos, pandemic scale | 5             | 1,000,000      | 5,000,000      | O(n), n~30    | ~50ms               |
| Wrong file, all mismatch  | 100,000       | 100,000        | 10,000,000,000 | O(n), n~30    | hours               |
| Wrong file, capped        | 20 (capped)   | 100,000        | 2,000,000      | O(n), n~30    | ~20ms               |

The "wrong file" scenario (all names unmatched) must be capped. When U exceeds a threshold, the error is systematic, not individual typos. The proposed cap: compute fuzzy suggestions for at most 20 unmatched names. Report all unmatched names in the batch error, but annotate only the first 20 with near-match suggestions.

Additional candidate filtering reduces M per query:

- Length filter: skip candidates differing in length by more than `max(3, len/3)` from the query. For typical taxon names, this eliminates 70-90% of candidates
- Score threshold: skip candidates scoring below 0.7 Jaro-Winkler similarity (too dissimilar to be useful)

With both filters on a 1M candidate set, effective M per query drops to 10K-100K.

#### Suggestion quality

String similarity does not understand domain structure. `A/California/07/2009` vs `A/California/08/2009` (Jaro-Winkler ~0.98) is a different strain, not a typo. Suggestions can surface misleading near-matches.

Mitigations:

- Suggestions are advisory only -- displayed in error messages, never used for matching
- Show the similarity score so the user can judge relevance
- Cap at top-3 per unmatched name to limit noise

### Axis 4: Name normalization

Unlike axes 1-3, this axis changes matching semantics: different inputs succeed that would otherwise fail, and different sequences attach to different nodes. Wrong normalization causes wrong data attachment, which produces wrong phylogenetic results silently.

Pros:

- Solves real-world mismatches that fuzzy suggestions can only diagnose but not fix (user still has to manually edit files)
- Ecosystem precedent: 2 of 9 surveyed tools normalize case (DendroPy, MrBayes), 1 normalizes characters (IQ-TREE), 2 normalize underscores (DendroPy, scikit-bio)
- Underscore-space equivalence in particular is the Newick standard, not an invention

Cons:

- Changes which data attaches to which node -- scientific correctness risk if normalization creates false matches
- Requires post-normalization collision detection (two distinct names folding to the same form)
- Adds CLI flags and configuration surface area
- Complicates reasoning about name identity throughout the pipeline (is `Sample_1` the same entity as `Sample 1`?)

Recommendation: 4a and 4b as opt-in CLI flags, off by default. Both require post-normalization collision detection. 4c deferred.

#### 4a: Case-insensitive matching (`--match-case-insensitive`)

Fold both tree and alignment names to lowercase before matching. Precedent: DendroPy (default, [[src](https://github.com/jeetsukumaran/DendroPy/blob/75bc2aea3450ef48a7aba5017551318f6586b0f6/src/dendropy/datamodel/taxonmodel.py#L538)]), MrBayes (always, [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/utils.c#L1707-L1726)]).

Pros:

- Solves the common `USA` vs `usa` class of mismatch
- Well-established precedent (2 tools do this by default)

Cons:

- Names that are intentionally case-distinct (`USA` and `usa` as separate taxa) would collide
- Must detect and error on post-folding collisions

#### 4b: Underscore-space equivalence (`--match-underscore-as-space`)

Treat `_` and ` ` as equivalent in name comparisons. Precedent: DendroPy (`preserve_underscores=False` default), scikit-bio (`convert_underscores=True` default), Newick standard (Felsenstein, 1986).

Two implementation approaches:

- At match time: normalize `_` to ` ` in both map keys and lookup keys before matching. Store original names for output. Scoped to the matching boundary -- name identity unchanged elsewhere in the pipeline
- At parse time: apply the Newick standard conversion in the Newick reader ([[src](https://github.com/neherlab/treetime/blob/ea8a21a74d8cbfd36466caa89f60e4736621de82/packages/util-newick/src/parse.rs)]). Makes `_`-to-space the default for all Newick inputs, matching DendroPy and scikit-bio. But changes name identity globally, not just at matching

Pros (match-time):

- Scoped -- only affects matching, not downstream name identity
- Explicit opt-in via CLI flag

Pros (parse-time):

- Newick-standard-compliant by default
- Consistent with DendroPy and scikit-bio behavior
- No CLI flag needed

Cons (match-time):

- Non-standard: Newick spec says underscores ARE spaces in unquoted labels; treating them as different is the deviation

Cons (parse-time):

- Changes name identity globally -- output files would contain spaces where input had underscores
- Breaks round-trip: `tree.nwk` -> parse -> write might produce different names
- Names with intentional underscores (not representing spaces) would be corrupted

Recommendation: match-time normalization as opt-in flag. The parse-time approach is what the Newick standard specifies, but the ecosystem consensus (7 of 9 tools keep underscores verbatim) suggests that users do not expect underscores to become spaces.

#### 4c: Character normalization (`--match-normalize-names`)

Replace non-alphanumeric characters (except `_`, `-`, `.`) with `_` in both sources before matching. Precedent: IQ-TREE `renameString()` ([[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/utils/tools.cpp#L503-L511)]).

Pros:

- Handles the broadest class of mismatches (`A/California/07/2009` vs `A_California_07_2009`)
- Useful for interop with IQ-TREE-generated trees

Cons:

- Highest false-match risk: many semantically distinct names collapse to the same normalized form
- Only 1 of 9 tools does this (IQ-TREE), and IQ-TREE applies it unconditionally (not opt-in)
- No concrete user request for this

Recommendation: defer until a concrete user need emerges.

## Generality across input sources

The proposal introduction mentions four input sources (FASTA, dates, traits, tree), but the analysis so far focuses on FASTA+tree. This section examines whether the matching problem is the same across all sources and whether a shared solution is appropriate.

### What varies across reconciliation sites

| Dimension                  | FASTA sequences                         | Date constraints                                         | Mugration traits                         |
| -------------------------- | --------------------------------------- | -------------------------------------------------------- | ---------------------------------------- |
| Tree scope                 | leaves only                             | all nodes (leaves + internal)                            | leaves only                              |
| External source type       | `Vec<FastaRecord>` (unindexed)          | `BTreeMap<String, Option<DateConstraint>>` (pre-indexed) | `BTreeMap<String, String>` (pre-indexed) |
| Missing leaf in external   | fill with ambiguous if <1/3, else error | mark `bad_branch`, warn                                  | hard error                               |
| Missing external in tree   | silently ignored                        | warn (unused constraints)                                | hard error                               |
| Duplicate risk in external | yes (`Vec` allows dupes)                | no (map keys unique)                                     | no (map keys unique)                     |
| Lookup method              | `.find()` linear scan                   | `BTreeMap::get` O(log n)                                 | `BTreeMap::get` O(log n)                 |

### What is the same

The core operation at every site:

1. Collect tree node names (from a configurable scope: leaves, all nodes)
2. Collect external data names (from a heterogeneous source)
3. Compute: matched pairs, unmatched-in-tree, unmatched-in-external, duplicates-in-tree, duplicates-in-external
4. Optionally compute fuzzy suggestions for unmatched names
5. Return the reconciliation result; let the caller decide policy (error, warn, fill)

This is a set-reconciliation operation parameterized by policy. The algorithm (index, diff, dedup, suggest) is identical; only the reaction to results differs.

### Current state: three independent reimplementations

- FASTA: worst implementation (linear scan, first-miss abort, no dedup)
- Dates: middle implementation (pre-indexed input, partial batch reporting via `warn_unused_date_constraints`, no bidirectional mismatch report)
- Mugration: best implementation (`IndexSet::difference` bidirectional batch, but no dedup because input is already a map, no fuzzy suggestions)

Each site reinvents the name-collection and diff logic. The dates site has its own `warn_unused_date_constraints` function that is essentially a one-directional set-difference. The mugration site's `validate_trait_names` is a bidirectional set-difference. Neither is reusable by the others.

### Architecture options

#### Option A: Fix each site independently

Apply the index+batch+dedup improvements to each site separately, duplicating the logic.

Pros:

- Each site can optimize for its specific types and constraints
- No abstraction overhead or indirection
- Simplest diff per site

Cons:

- Three copies of the same algorithm with the same bugs and the same future improvements
- New reconciliation sites (config file partitions, multi-segment genomes) would need a fourth copy
- Normalization and fuzzy logic duplicated or absent from some sites

#### Option B: Shared reconciliation function, type-erased

Extract the algorithm into a shared function that operates on two sets of `&str` names. The function knows nothing about `FastaRecord`, `DateConstraint`, or trait values -- it only reconciles name sets and returns the result. Each caller maps its domain types to name sets before calling, and interprets the result (error/warn/fill) after.

```rust
pub struct ReconciliationResult<'a> {
    pub matched_tree_names: Vec<&'a str>,
    pub unmatched_in_tree: Vec<&'a str>,
    pub unmatched_in_external: Vec<&'a str>,
    pub duplicate_tree_names: Vec<(&'a str, usize)>,
    pub duplicate_external_names: Vec<(&'a str, usize)>,
    pub suggestions: BTreeMap<&'a str, Vec<(&'a str, f64)>>,
}

pub struct ReconciliationConfig {
    pub case_insensitive: bool,
    pub underscore_as_space: bool,
    pub max_suggestions_per_name: usize,
    pub max_names_with_suggestions: usize,
    pub suggestion_threshold: f64,
}

pub fn reconcile_names<'a>(
    tree_names: &[&'a str],
    external_names: &[&'a str],
    config: &ReconciliationConfig,
) -> ReconciliationResult<'a>;
```

Pros:

- One implementation of index, diff, dedup, fuzzy -- tested and maintained once
- All sites automatically benefit from normalization options and fuzzy suggestions
- New reconciliation sites (config file, multi-segment) call the same function
- Type-erased: no dependency on `FastaRecord` or any domain type -- lives in `treetime-utils`
- Policy stays at the call site where domain knowledge lives (dates: warn+mark-bad; FASTA: fill-or-error; mugration: hard error)

Cons:

- Extra name-extraction step at each call site (collect `&str` names from domain types before calling)
- Dates site operates on all nodes, not just leaves -- caller must filter the tree scope before calling
- For pre-indexed sources (`BTreeMap` dates, traits), the function rebuilds an internal `HashMap` from the provided `&[&str]`, which is redundant. Minor cost (O(n) on already-indexed data)

#### Option C: Trait-based reconciliation

Define a `NameSource` trait that abstracts over different input types:

```rust
pub trait NameSource {
    fn names(&self) -> Vec<&str>;
    fn has_duplicates(&self) -> bool;
}
```

Implement for `Vec<FastaRecord>`, `BTreeMap<String, T>`, etc. The reconciliation function accepts `impl NameSource` for both sides.

Pros:

- Type-safe: each source type implements its own name extraction
- Can optimize: `BTreeMap` impl knows keys are unique, skips duplicate check

Cons:

- Trait abstraction for a function that runs once at startup -- overengineered for the use case
- The trait interface is trivial (just `.names()`) -- a closure or `.iter().map()` achieves the same thing without a trait definition
- Harder to test: need to construct trait objects instead of plain `&[&str]`

### Recommendation: Option B

The reconciliation algorithm is a pure function on two name sets. The domain-specific types (`FastaRecord`, `DateConstraint`, trait `String`) are irrelevant to the matching logic -- only the names matter. Option B separates the algorithm from the policy cleanly:

- The shared function in `treetime-utils` handles: indexing, dedup detection, bidirectional set-difference, fuzzy suggestions, normalization
- Each call site handles: extracting names from domain types, choosing tree scope (leaves vs all), interpreting the result (error vs warn vs fill)

This replaces three ad-hoc reimplementations with one tested function, and naturally extends to future input sources.

Proposed location: `packages/treetime-utils/src/reconcile.rs`.

### How input formats affect reconciliation

Not all input paths require name matching. Unified formats embed data in the tree structure, bypassing the reconciliation step entirely. This table maps each format to whether it needs reconciliation for each data type:

| Input format           | Sequences                                        | Dates                       | Traits                          | Name matching needed?                                                                                                         |
| ---------------------- | ------------------------------------------------ | --------------------------- | ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| Newick + FASTA + TSV   | separate file                                    | separate TSV                | separate TSV                    | yes -- all three require name reconciliation                                                                                  |
| Auspice JSON           | mutations on branches (`branch_attrs.mutations`) | `num_date` in node attrs    | `country` etc. in node attrs    | no -- all data embedded in tree nodes/edges                                                                                   |
| UShER MAT protobuf     | mutations in preorder                            | not carried                 | not carried                     | no for sequences; yes if dates/traits needed (separate TSV)                                                                   |
| PhyloXML               | per-clade `<sequence>` elements                  | per-clade `<date>` elements | per-clade `<property>` elements | no -- all data embedded in tree                                                                                               |
| Nexus + FASTA + TSV    | tree only (TRANSLATE resolves names internally)  | separate TSV                | separate TSV                    | yes for sequences/dates/traits, but Nexus TRANSLATE table handles tree-internal name aliasing                                 |
| Config file (proposed) | per-partition alignment paths                    | per config                  | per config                      | depends on whether config includes explicit name mappings (see [config-file-multi-partition](config-file-multi-partition.md)) |

Key observations:

- Auspice JSON and PhyloXML are fully self-contained: sequences, dates, and traits are all node/edge attributes. No name matching needed. The reader populates the graph directly
- UShER MAT carries only mutations. If a timetree analysis needs dates, a separate TSV is still required, and name matching applies to the dates-to-tree reconciliation even though sequences are embedded
- The Nexus TRANSLATE table is a name-mapping mechanism built into the format: numeric tokens in the Newick string map to full taxon names. The v1 Nexus parser (`util-newick/src/nexus.rs`) already resolves these during parsing. This is not the same as the proposal's normalization (axis 4) -- it's format-level aliasing, handled before matching begins
- The proposed config file format could include explicit name mappings, making reconciliation configurable per-partition

The reconciliation function from Option B fits naturally: it runs only when the input path requires name matching (Newick+FASTA, or UShER MAT + dates TSV). Unified-format readers skip it entirely. This makes the function a component of the Newick+separate-files input path, not a universal pipeline stage.

## Validation plan

### Correctness

- Existing tests with matching names pass unchanged
- Duplicate leaf names in tree: error listing all duplicates
- Duplicate FASTA names: error listing all duplicates
- Single mismatch: error with fuzzy suggestions
- Multiple mismatches: batch error listing all, suggestions for first 20
- Case mismatch (`USA` vs `usa`): exact match fails; suggestion shows the case-different name; `--match-case-insensitive` makes it succeed
- Underscore/space swap (`Sample_1` vs `Sample 1`): exact match fails; suggestion shows the variant; `--match-underscore-as-space` makes it succeed
- Post-normalization collision: error even with normalization flags enabled

### Performance targets

- 1K sequences: indexing + matching <1ms
- 100K sequences: indexing + matching <100ms
- 1M sequences: indexing + matching <1s
- 10 unmatched among 100K candidates with fuzzy suggestions: <50ms
- 10 unmatched among 1M candidates with fuzzy suggestions: <500ms
- 100K unmatched (wrong file): capped at 20 fuzzy, rest reported without suggestions
- Benchmark with `criterion` at 1K, 10K, 100K, 1M name sets

### Regression

- Golden master tests for ancestral, timetree, optimize, mugration pass unchanged
- Error message format changes: update tests that check exact error text

## Implementation order

1. Axis 1 + 2b (index + duplicate detection): combined, since index construction naturally detects duplicates. Fixes the quadratic performance regression and the silent-wrong-data defect
2. Axis 2a (batch mismatch reporting): requires the index from step 1. Brings sequence attachment to parity with the mugration validator internally and IQ-TREE's `setAlignment()` externally
3. Axis 3 (fuzzy suggestions): builds on the batch mismatch data from step 2
4. Axis 4a-4b (normalization): independent CLI flags, can be done in any order after step 1

Steps 1-2 are correctness fixes. Step 3 is a diagnostic improvement. Step 4 changes matching semantics.

## Related

### Issues

- [kb/issues/M-io-sequence-attachment-quadratic.md](../issues/M-io-sequence-attachment-quadratic.md) -- O(n^2) performance, resolved by Axis 1
- [kb/issues/M-io-sequence-name-matching-unreliable.md](../issues/M-io-sequence-name-matching-unreliable.md) -- source issue for this proposal
- [kb/issues/M-dates-column-auto-detection-gaps.md](../issues/M-dates-column-auto-detection-gaps.md) -- similar name reconciliation defects in CSV metadata readers
- [kb/issues/N-io-large-dataset-memory-constraint.md](../issues/N-io-large-dataset-memory-constraint.md) -- large-dataset memory considerations; the name index adds negligible overhead relative to sequence data

### Tickets

- [kb/tickets/io-sequence-name-matching-unreliable.md](../tickets/io-sequence-name-matching-unreliable.md) -- implementation ticket derived from this proposal

### Proposals

- [kb/proposals/config-file-multi-partition.md](config-file-multi-partition.md) -- configuration file format that could include explicit name mappings as an alternative to fuzzy matching
- [kb/proposals/unified-input-format-support.md](unified-input-format-support.md) -- unified input formats (Auspice JSON, UShER MAT) bypass name matching entirely by embedding sequence data in the tree structure

### Decisions

- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md) -- documents the Newick underscore-space standard and annotation dialect landscape
