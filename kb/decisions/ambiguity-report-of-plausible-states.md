# The ambiguity report names the states near the most likely state

The optional ambiguity report of `ancestral`, `timetree`, and `mugration` names, for each reconstructed node and position, every state whose posterior probability is at most a threshold below the largest probability. The threshold is the value of `--ambiguity-threshold`; without it the report is off. The outputs that carry the report are selected with `--ambiguity-outputs`, by default every output of the run that can carry it. v0 and augur write no such report.

**Type**: Scientific output; new in v1.

**Status**: Approved on 2026-10-21.

## Rule

- The plausible states of a posterior $p$ with threshold $t$ are the states $k$ with $\max_j p_j - p_k \le t + 10^{-10}$, the tie tolerance of [kb/decisions/posterior-most-likely-state-tie-rule.md](posterior-most-likely-state-tie-rule.md). The set therefore always contains the most likely state, and with $t = 0$ it names the tied states only. `AmbiguityThreshold::plausible_states()` in [packages/treetime/src/partition/marginal/posterior_states.rs](../../packages/treetime/src/partition/marginal/posterior_states.rs) owns the rule for traits and sequences
- The threshold is a difference of probabilities from 0 to 1. There is no default value. RAxML-NG uses $0.5/n$ for $n$ states, for example 0.125 for nucleotides
- A node with more than one plausible state is reported. The most likely state, the mutations, the branch labels, and every other output stay as without the report

| Output | Carrier name | Report |
| --- | --- | --- |
| Reconstructed nucleotide FASTA (`ancestral`, `timetree`) | `reconstructed-nuc-fasta` | IUPAC code of the plausible states at the position, `N` for all four |
| Reconstructed amino-acid FASTA (`ancestral`) | `reconstructed-aa-fasta` | `B` for `D` and `N`, `Z` for `E` and `Q`, `J` for `I` and `L`, otherwise `X` |
| Node sequences of the augur node data (`ancestral`) | `augur-node-data` | the same codes as the nucleotide FASTA |
| Traits CSV (`mugration`) | `traits-csv` | column `<attribute>_ambiguous`: the plausible state names joined by `\|` in the order of the state names, empty for an unambiguous node |
| Augur node data (`mugration`) | `augur-node-data` | node field `<attribute>_ambiguous` with the same text, only for ambiguous nodes |
| Auspice JSON (`mugration`) | `auspice` | node attribute `<attribute>_ambiguous` with the same text as `value`, only for ambiguous nodes. Auspice shows every node attribute in the panel of a clicked node |

- Sequences: internal nodes have a report at every position outside gaps and unknown regions, leaves only at the positions that `--impute-missing-data` imputes. A node whose sequence `--sample-from-profile` draws has no report, and `--method-anc=parsimony` gives no report, because Fitch parsimony computes no posterior. The sparse representation reads the profile of each fixed state for the positions that carry that state. Both representations give the same report, except at the unknown columns of an imputed tip below a variable column of its parent, where the sparse backend has no posterior of the tip ([kb/issues/M-ancestral-sparse-imputation-of-unknown-tip-columns-ignores-parent-posterior.md](../issues/M-ancestral-sparse-imputation-of-unknown-tip-columns-ignores-parent-posterior.md))
- A carrier that the run does not write gives a warning when `--ambiguity-outputs` names it

## Reference behavior

RAxML-NG writes ambiguous ancestral states in `AncestralStates::ml_state()` (`.repos/raxml-ng/src/AncestralStates.cpp` lines 60-95 at revision `d396351ee1263b704adb135edd6a4a84522dbd5b`): a state is in the set when $p_k + 0.5/n > \max_j p_j$, and the set is written as the code of the set when the alphabet has one, otherwise as the most likely state. v1 differs in three points: the threshold is a value of the run, the comparison includes the tie tolerance, and a set without a code is written as `N` or `X`, so the report never shows a state as certain when it is not. RAxML-NG has no trait reconstruction.

## Why

A single most likely state hides that the data cannot separate two or more states. The report shows this in every output whose format can hold it, as codes in sequences and as extra fields in tables and JSON documents, and leaves the most likely state, the mutations, and the formats of the other outputs unchanged.
