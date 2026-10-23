# Input name matching and validation for tree-data reconciliation

Analysis commands pair entries of separate inputs with tree nodes by name: FASTA records, metadata dates, and trait values. The pairing is exact and case-sensitive, so a name that differs from its leaf only in case, whitespace, or punctuation leaves the leaf without data. This proposal analyzes two open design axes of that pairing: suggestions for unmatched names, and opt-in name normalization. An ecosystem survey of nine phylogenetics tools informs both.

## Current state

- **One pairing function**: `fn pair_by_name()` [packages/treetime-graph/src/pair_by_name.rs#L8-L49](../../packages/treetime-graph/src/pair_by_name.rs#L8-L49) pairs entries with candidate nodes for every input. It builds one index of entry names while it reads the entries, looks the input name of each node up in that index, and returns the entries that match no node. The callers are `fn pair_alignment()` for FASTA records [packages/app-commands/src/commands/shared/alignment.rs#L73](../../packages/app-commands/src/commands/shared/alignment.rs#L73), `fn read_input_dates()` for dates [packages/app-commands/src/commands/shared/dates_input.rs#L17](../../packages/app-commands/src/commands/shared/dates_input.rs#L17), and `fn pair_trait_rows()` for trait values [packages/treetime/src/mugration/trait_pairing.rs#L7](../../packages/treetime/src/mugration/trait_pairing.rs#L7)
- **Repeated names**: the tree reader gives each later node of a repeated name the unique name `<name>_<k>`, and the pairing uses the input names: the first entry of a name wins, and every node with that input name receives it. Typed run warnings name the repeated names. [kb/decisions/duplicate-names-warned-ids-from-input-order.md](../decisions/duplicate-names-warned-ids-from-input-order.md) records this rule
- **Mismatch reports**: a leaf without a sequence or a trait value gets one warning each and a count warning; dates count such leaves in the date summary. Dates and traits warn about entries without a node with up to 10 names; FASTA records without a leaf get no message. [kb/issues/N-io-missing-name-policy-differs-across-subsystems.md](../issues/N-io-missing-name-policy-differs-across-subsystems.md) tracks these differences, with an open decision
- **No suggestions and no normalization**: the messages name the unmatched names but no near-match. The only fuzzy matching in v1 is `fn suggestion_suffix()` [packages/app-commands/src/config/suggest.rs#L4](../../packages/app-commands/src/config/suggest.rs#L4), which suggests config keys and values by Levenshtein distance through `strsim`. [kb/issues/N-io-sequence-name-matching-unreliable.md](../issues/N-io-sequence-name-matching-unreliable.md) tracks the user-facing defects

Inputs that carry their data inside the tree (Auspice JSON, PhyloXML) need no name pairing. UShER MAT carries mutations only, so dates and traits for it still pair by name. See [unified-input-format-support.md](unified-input-format-support.md).

## Ecosystem survey

Nine phylogenetics tools were inspected at the source level to learn how they reconcile tree names with alignment names. For each tool, the comparison function, the duplicate detection, and the error reporting were traced in source. Every claim below links to the source location.

### Matching semantics

All tools use exact string equality for the primary lookup. Where normalization exists, it is applied symmetrically to both sides as a preprocessing step before the exact match -- not as a fuzzy fallback. This distinction matters: normalization changes which data attaches to which tree node and therefore affects scientific results. Suggestions for unmatched names (axis 1 of this proposal) are diagnostic only and do not change matching behavior.

| Tool        | Comparison                                             | Normalization                             | Source                                                                                                                                                                                                                                                                                                 |
| ----------- | ------------------------------------------------------ | ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| v0 TreeTime | `name not in dict` (exact)                             | none                                      | [[src](../../packages/legacy/treetime/treetime/treeanc.py#L395-L444)]                                                                                                                                                     |
| augur       | delegates to v0                                        | space-split for IQ-TREE compat only       | [[src](https://github.com/nextstrain/augur/blob/024292af6daf1f8aae7e8d690bd2f944e1c7322d/augur/align.py#L197-L211)]                                                                                                                                                                                    |
| IQ-TREE 2   | exact after `renameString()`                           | non-alphanum -> `_` on both sides         | [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/utils/tools.cpp#L503-L511)], [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/tree/phylotree.cpp#L349-L383)]                                                                |
| BEAST2      | `String.equals()` (exact)                              | none                                      | [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/alignment/Alignment.java#L321)], [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/TreeParser.java#L424)]         |
| BEAST v1    | `String.equals()` (exact)                              | FASTA header trim only                    | [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/util/TaxonList.java#L114)], [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evomodel/treelikelihood/BeagleTreeLikelihood.java#L490)] |
| MrBayes     | exact after `tolower()`                                | case fold on both sides (always)          | [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/utils.c#L1707-L1726)]                                                                                                                                                                                    |
| DendroPy    | `str.lower()` when `is_case_sensitive=False` (default) | unquoted `_` -> space per Newick standard | [[src](https://github.com/jeetsukumaran/DendroPy/blob/75bc2aea3450ef48a7aba5017551318f6586b0f6/src/dendropy/datamodel/taxonmodel.py#L538)]                                                                                                                                                             |
| gotree      | Go map lookup (exact)                                  | none                                      | [[src](https://github.com/evolbioinfo/gotree/blob/83cca67fdcf64af47f6f3318d3eebdcb417d45a4/tree/tree.go#L454-L468)]                                                                                                                                                                                    |
| ETE 4       | Python dict lookup (exact)                             | none                                      | [[src](https://github.com/etetoolkit/ete/blob/9562dfb6a02795dfda3975be5025f83b8dc884b1/ete4/phylo/phylotree.py#L346-L369)]                                                                                                                                                                             |

IQ-TREE's `renameString()` is the most aggressive normalizer: it replaces any character outside `[A-Za-z0-9_.-]` with `_`, applied to both tree and alignment names at parse time. MrBayes applies `tolower()` unconditionally through its central `StrCmpCaseInsensitive()` function used in all name lookups. DendroPy is the only tool implementing the Newick-standard underscore-to-space rule (`preserve_underscores=False` default in its Newick reader).

### Duplicate detection

| Tool        | Tree duplicates                   | Alignment duplicates                 | Source                                                                                                                                                                                                                                                                                              |
| ----------- | --------------------------------- | ------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MrBayes     | hard error (all contexts)         | hard error                           | [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/command.c#L7234-L7285)]                                                                                                                                                                               |
| gotree      | hard error on index build         | n/a                                  | [[src](https://github.com/evolbioinfo/gotree/blob/83cca67fdcf64af47f6f3318d3eebdcb417d45a4/tree/tree.go#L454-L468)]                                                                                                                                                                                 |
| BEAST2      | hard error (first)                | hard error (first)                   | [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/TreeParser.java#L462-L471)], [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/alignment/Alignment.java#L321)] |
| IQ-TREE 2   | none                              | batch (sorted adjacent scan)         | [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/alignment/alignment.cpp#L136-L161)]                                                                                                                                                                           |
| augur       | none                              | batch (sorted set, all dupes listed) | [[src](https://github.com/nextstrain/augur/blob/024292af6daf1f8aae7e8d690bd2f944e1c7322d/augur/io/sequences.py#L506-L540)]                                                                                                                                                                          |
| BEAST v1    | silent dedup in `Taxa.addTaxon()` | none                                 | [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/util/TaxonList.java#L114)]                                                                                                                                                            |
| v0 TreeTime | none                              | silent (dict key overwrite)          | [[src](../../packages/legacy/treetime/treetime/sequence_data.py#L180-L194)]                                                                                                                                                       |
| ETE 4       | error in EvolTree only            | silent auto-rename (`2_name`)        | [[src](https://github.com/etetoolkit/ete/blob/9562dfb6a02795dfda3975be5025f83b8dc884b1/ete4/parser/fasta.py#L8-L50)]                                                                                                                                                                                |

MrBayes and gotree are the strictest: any duplicate name is a hard error. BEAST2 detects duplicates in both sources but reports only the first. IQ-TREE and augur detect alignment duplicates in batch but not tree duplicates. v0, BEAST v1, and ETE handle duplicates silently (overwrite, dedup, or rename), masking data errors.

### Batch error reporting

| Tool        | Behavior                                                                                         | Source                                                                                                                                             |
| ----------- | ------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| IQ-TREE 2   | bidirectional batch: all alignment-not-in-tree + all tree-not-in-alignment reported before abort | [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/tree/phylotree.cpp#L349-L383)]                               |
| augur       | batch for FASTA duplicates (sorted list); one-at-a-time elsewhere                                | [[src](https://github.com/nextstrain/augur/blob/024292af6daf1f8aae7e8d690bd2f944e1c7322d/augur/io/sequences.py#L506-L540)]                         |
| BEAST2      | set-difference warnings (non-fatal, both directions); one-at-a-time for hard errors              | [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/TreeParser.java#L162-L181)]  |
| v0 TreeTime | per-leaf warning + cumulative count + abort when >1/3 missing                                    | [[src](../../packages/legacy/treetime/treetime/treeanc.py#L395-L444)] |
| All others  | one-at-a-time (abort on first mismatch)                                                          | see matching semantics table                                                                                                                       |

IQ-TREE's `setAlignment()` is the best implementation: it iterates both directions (alignment -> tree, tree -> alignment), collects all mismatches, and reports them in a single error. This is the external model for batch reporting.

### Fuzzy suggestions

No surveyed tool provides fuzzy matching suggestions for unmatched names. This would be novel in the phylogenetics ecosystem.

### Newick underscore-space standard

The Newick specification (Felsenstein, 1986) states that underscores in unquoted labels represent spaces: "it is assumed that an underscore character stands for a blank; any of these in a name will be converted to a blank when it is read in" ([source](https://phylipweb.github.io/phylip/newicktree.html)). See also [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md).

DendroPy and scikit-bio implement this conversion by default. All other surveyed tools (v0, augur, BEAST, IQ-TREE, gotree, ETE, MrBayes) keep underscores verbatim, and so does v1: the Newick reader of `util-newick` has the read option `underscores_as_spaces` ([packages/util-newick/src/read/options.rs#L11](../../packages/util-newick/src/read/options.rs#L11)), and TreeTime sets it to `false` ([packages/treetime-io/src/nwk.rs#L38](../../packages/treetime-io/src/nwk.rs#L38)).

## Design axes

The two axes are independent. Each can be adopted, deferred, or rejected separately.

### Axis 1: Suggestions for unmatched names

For each unmatched name, compute the string similarity to the names of the other input and show the closest names in the warning. Suggestions are advisory only: matching stays exact. No surveyed tool provides this.

Example warning:

```
No sequence found for leaf 'A/California/07/2009' in alignment 'aln.fasta'; treating it as fully ambiguous (missing data). Closest sequence names: 'A_California_07_2009' (0.93)
```

Pros:

- A message that names the closest names shows the cause (case, underscore, truncation) without a manual diff of the inputs
- No risk to results, because suggestions never change the pairing
- No new dependency: `strsim` is already a direct dependency of `app-commands`, and the config suggestions use it

Cons:

- The search is O(U * M) for U unmatched names and M candidates, so it needs a cap for large inputs (see "Performance at scale")
- Names that are close but distinct produce misleading suggestions: `A/California/07/2009` and `A/California/08/2009` are different strains, not a typo
- Message formatting becomes more complex

#### Algorithm

Taxon names are short strings (10 to 100 characters, mostly ASCII). Candidate metrics:

| Metric              | Cost per pair | Properties                                                                                                                                       |
| ------------------- | ------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| Jaro-Winkler        | O(n)          | Weights the shared prefix. Phylogenetic names share long prefixes (`A/California/07/2009` vs `A/California/7/2009`), so this metric fits them |
| Levenshtein         | O(n\*m)       | Insertions, deletions, substitutions. Easy to interpret. The config suggestions use it                                                           |
| Damerau-Levenshtein | O(n\*m)       | Adds transpositions. Finds `AB` vs `BA` swaps                                                                                                    |

Jaro-Winkler costs O(n) per pair, gives a score from 0 to 1, and its prefix weight fits names whose variable part is at the end. Levenshtein keeps one metric for all suggestions of the CLI. `strsim` provides both. For more throughput, `rapidfuzz` (batch comparison with a score cutoff), `editdistancek` (bounded Levenshtein), and `triple_accel` (SIMD) exist.

#### Performance at scale

The search runs only for unmatched names. Its cost is U * M * c, with U unmatched names, M candidates, and c the cost per pair.

| Scenario                  | U (unmatched) | M (candidates) | Pairs          | Cost per pair | Estimated wall time |
| ------------------------- | ------------- | -------------- | -------------- | ------------- | ------------------- |
| Few typos, small dataset  | 5             | 1,000          | 5,000          | O(n), n~30    | <1ms                |
| Few typos, large dataset  | 5             | 100,000        | 500,000        | O(n), n~30    | ~5ms                |
| Few typos, pandemic scale | 5             | 1,000,000      | 5,000,000      | O(n), n~30    | ~50ms               |
| Wrong file, all mismatch  | 100,000       | 100,000        | 10,000,000,000 | O(n), n~30    | hours               |
| Wrong file, capped        | 20 (capped)   | 100,000        | 2,000,000      | O(n), n~30    | ~20ms               |

When all names are unmatched, the cause is systematic (a wrong file), and suggestions do not help. The proposed cap: suggestions for at most 20 unmatched names. Two filters reduce M per name:

- Length filter: skip candidates whose length differs by more than `max(3, len/3)`, which removes 70% to 90% of the candidates for typical taxon names
- Score threshold: skip candidates below a Jaro-Winkler similarity of 0.7

#### Suggestion quality

Mitigations for misleading suggestions:

- Suggestions appear in messages only and never change the pairing
- Each suggestion shows its similarity score, so the user can judge it
- At most 3 suggestions per unmatched name

### Axis 2: Name normalization

This axis changes which entries pair with which nodes: inputs that fail without it succeed with it. A wrong normalization attaches wrong data to a node and changes results without a message.

Pros:

- Solves mismatches that suggestions can only show; without normalization the user must edit the input files
- Ecosystem precedent: 2 of 9 surveyed tools fold case (DendroPy, MrBayes), 1 replaces characters (IQ-TREE), and the Newick standard treats underscores as spaces (DendroPy and scikit-bio follow it)

Cons:

- A normalization that joins distinct names attaches data to the wrong node
- Folding can make distinct names equal. Under [kb/decisions/duplicate-names-warned-ids-from-input-order.md](../decisions/duplicate-names-warned-ids-from-input-order.md) such names are repeated names: the first entry wins and a warning names them. Stopping the run instead needs an exception to that rule
- Adds CLI flags and configuration
- Makes name identity harder to reason about (is `Sample_1` the same sample as `Sample 1`?)

#### 2a: Case-insensitive matching (`--match-case-insensitive`)

Fold tree names and entry names to lowercase before the pairing. Precedent: DendroPy (default, [[src](https://github.com/jeetsukumaran/DendroPy/blob/75bc2aea3450ef48a7aba5017551318f6586b0f6/src/dendropy/datamodel/taxonmodel.py#L538)]), MrBayes (always, [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/utils.c#L1707-L1726)]).

- Pro: solves the frequent `USA` vs `usa` mismatch, with precedent in 2 tools
- Con: names that differ only in case on purpose become one name

#### 2b: Underscore-space equivalence (`--match-underscore-as-space`)

Treat `_` and space as equal in the pairing. Precedent: DendroPy (`preserve_underscores=False` by default), scikit-bio (`convert_underscores=True` by default), and the Newick standard (Felsenstein, 1986).

Two places can apply it:

- At pairing time: replace `_` with a space in the index keys and in the looked-up names, and keep the original names for the outputs. Only the pairing changes; names stay the same everywhere else
- At read time: turn on `underscores_as_spaces` of the Newick reader. This follows the Newick standard, but it changes the names in every output (spaces where the input had underscores), breaks the round trip of names through TreeTime, and corrupts names whose underscores are not spaces

Recommendation: pairing time, as an opt-in flag. The read-time conversion is what the Newick standard specifies, but 7 of 9 surveyed tools keep underscores, so users do not expect them to become spaces.

#### 2c: Character normalization (`--match-normalize-names`)

Replace every character other than letters, digits, `_`, `-`, and `.` with `_` in both inputs before the pairing. Precedent: IQ-TREE `renameString()` ([[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/utils/tools.cpp#L503-L511)]).

- Pro: covers the widest class of mismatches (`A/California/07/2009` vs `A_California_07_2009`), and helps with trees written by IQ-TREE
- Con: the highest risk of joining distinct names; only IQ-TREE does this, and it does it always

Recommendation: 2a and 2b as opt-in flags, off by default. Defer 2c until a user needs it.

## Validation plan

- Existing tests with matching names pass unchanged, and the golden masters of `ancestral`, `timetree`, `optimize`, and `mugration` pass unchanged
- One unmatched name: the warning shows the closest names with scores
- Many unmatched names: suggestions for the first 20 only, and every unmatched name still in the report
- Case mismatch (`USA` vs `usa`): the exact pairing fails and the suggestion shows the other case; `--match-case-insensitive` pairs them
- Underscore and space (`Sample_1` vs `Sample 1`): the exact pairing fails and the suggestion shows the variant; `--match-underscore-as-space` pairs them, and the outputs keep the input names
- Names that folding makes equal: handled as repeated names, with a warning
- Performance, measured with `criterion`: 10 unmatched names among 100,000 candidates in under 50 ms, and among 1,000,000 candidates in under 500 ms

## Implementation order

1. Axis 1 (suggestions): changes messages only
2. Axis 2a and 2b (normalization flags): independent of each other; they change which entries pair with which nodes

## Related

### Issues

- [kb/issues/N-io-sequence-name-matching-unreliable.md](../issues/N-io-sequence-name-matching-unreliable.md): user-facing name-matching defects, the source issue of this proposal
- [kb/issues/N-io-missing-name-policy-differs-across-subsystems.md](../issues/N-io-missing-name-policy-differs-across-subsystems.md): what each input does with unmatched names
- [kb/issues/M-mugration-column-detection-no-positional-fallback.md](../issues/M-mugration-column-detection-no-positional-fallback.md): mugration column detection lacks the positional fallback of v0

### Proposals

- [config-file-multi-partition.md](config-file-multi-partition.md): a configuration file that could hold explicit name mappings
- [unified-input-format-support.md](unified-input-format-support.md): input formats that carry their data inside the tree and need no name pairing

### Decisions

- [kb/decisions/duplicate-names-warned-ids-from-input-order.md](../decisions/duplicate-names-warned-ids-from-input-order.md): repeated names are warned about, and identity comes from input order
- [kb/decisions/multi-format-tree-io.md](../decisions/multi-format-tree-io.md): the Newick underscore-space standard and the annotation dialects
