# Partition Configuration for TreeTime v1

## 1. Current state

Every command (ancestral, optimize, timetree, prune) creates exactly one sequence partition inline with hardcoded partition index `0`. The creation sequence is procedural code repeated per command: `create_fitch_partition` -> GTR resolution -> `log_gtr` -> conversion to marginal dense or sparse. Dense vs sparse is decided by a `--dense` CLI flag; when omitted, `infer_dense()` is called, which always returns `false` (stub).

The multi-partition infrastructure exists and works. Mugration adds a discrete trait partition alongside the sequence partition, demonstrating that the `Vec<Arc<RwLock<dyn PartitionTimetreeAll>>>` architecture handles heterogeneous partition types. But no command creates more than one sequence partition. No `PartitionConfig` type exists. Partition creation is not separated from command orchestration.

v0 has no explicit partition concept. `TreeAnc` owns a single `self.gtr` and single `self.data`. One alignment, one model.

## 2. Configuration methods

Five distinct methods for specifying what partitions to create have been discussed in the project. They are not mutually exclusive: a partition configuration layer would accept input from any of these sources and produce a uniform `Vec<PartitionSpec>`.

### 2.1 `--dense` flag (current default, v0 parity)

Single `--dense: Option<bool>` CLI flag on ancestral, optimize, and timetree commands. When `Some(true)`, creates one dense partition. When `Some(false)` or `None`, creates one sparse partition. Prune hardcodes sparse with JC69.

This is v0 parity: one alignment, one model, one partition.

**Sources:** `commands/ancestral/args.rs`, `commands/optimize/args.rs`, `commands/timetree/args.rs`

### 2.2 Auto-detection heuristic (`infer_dense`)

When `--dense` is omitted, `infer_dense()` in `partition/algo/infer_dense.rs` determines dense vs sparse. The stub always returns `false`. The planned heuristic: select dense when tree branches are long enough that most positions become variable (expected mutations per branch approaches sequence length). Below that threshold, sparse saves memory and compute.

**Sources:** `kb/issues/N-representation-infer-dense-stub.md`, `kb/tickets/representation-implement-infer-dense-heuristic.md`

### 2.3 Multiple alignment files

The design document `_raw/optimize.md` specifies "alignment(s) (corresponding to partitions)" and notes that "for more complex inputs (multiple alignments and models along with discrete characters) a config file format may be needed." Currently only single `--aln` is accepted. Multiple alignments would create one partition per file.

Ambiguity: when two alignment files are given, are they different genes on the same tree (multi-partition) or replicate alignments? The config file approach (S2.5) resolves this.

**Sources:** `kb/issues/N-optimize-multi-alignment-input.md`, `_raw/optimize.md`

### 2.4 `--segment` flag (multi-segment genomes)

For segmented viruses (influenza: 8 segments, same alphabet, independent evolution, shared tree), a `--segment` flag would load segment-aware FASTA and create one partition per segment. Branch `worktree/feat/multi-segment-genome-input` has 9 commits implementing this, unmerged. The design document `_raw/sequence_evolution.md` asks: "For flu, genomes come in segments -- should these sequences be saved as list of multiple sequences, or concatenated?"

v0 handles segments via concatenation.

**Sources:** `kb/issues/N-io-multi-segment-genome-input.md`, `_raw/sequence_evolution.md`

### 2.5 Config file (YAML/JSON)

A full proposal exists for a structured config file specifying:

- Partition list: each entry has alignment path, alphabet, model name or model file, dense/sparse mode
- Shared parameters: tree path, clock model, coalescent settings
- Output specification

This handles all use cases: multi-segment genomes, mixed data types (nucleotide + discrete location), per-partition model assignment, model file references for two-stage estimation. The config file complements CLI flags (simple cases use CLI, complex cases use config).

Prior art surveyed in the proposal:

| Tool     | Format                                   | Partition specification                         |
| -------- | ---------------------------------------- | ----------------------------------------------- |
| NEXUS    | SETS block                               | Character ranges within concatenated alignment  |
| IQ-TREE  | `-spp`/`-q`/`-sp` flags + partition file | Column ranges, three edge-linking modes         |
| RAxML-NG | Partition file                           | `MODEL, name = start-end` per line              |
| BEAST    | XML                                      | Full model hierarchy, link/unlink per partition |

Open design questions: YAML vs JSON vs TOML; config replaces or complements CLI; edge-linked vs edge-unlinked partitions.

Existing unmerged branches: `worktree/feat/optimize-multi-alignment` (config parsing + multi-partition tests).

**Sources:** `kb/proposals/config-file-multi-partition.md`

### 2.6 Unified input formats

Auspice JSON and UShER MAT files contain embedded sequence data. Loading these formats could populate partitions directly without separate alignment files.

**Sources:** `kb/proposals/unified-input-format-support.md`

## 3. Automatic partitioning (heuristic derivation)

Beyond user specification, partitions can be derived automatically from alignment properties. These methods take one alignment and produce a `position -> partition_index` map.

### 3.1 Rate-based partitioning (mutation count quantiles)

Run Fitch parsimony on the full alignment. Count mutations per site. Sort sites by mutation count. Split into $q$ quantile bins. Sites in each bin form one partition. Related to k-means site-rate clustering (Frandsen et al. 2015).

This creates partitions where within-partition rate variation is minimized, improving GTR model fit per partition.

**Sources:** `kb/reports/auto-partitioning.md` S7.1

### 3.2 Codon-position partitioning (from GFF annotation)

Parse GFF annotations to identify CDS regions. Split protein-coding positions by codon position (1st, 2nd, 3rd). Non-coding positions form a 4th partition. Can further isolate specific genes (e.g., Spike protein) into their own 1st/2nd/3rd codon-position partitions.

Biologically motivated: 3rd codon positions evolve faster (wobble), 1st/2nd positions are more constrained. Mixing them in one GTR model distorts both rate and frequency estimates.

Feature inventory status: `[ ]` not implemented.

**Sources:** `kb/reports/auto-partitioning.md` S7.1, `kb/features/representation.md`

## 4. Model configuration (how models relate to partitions)

### 4.1 GTR family grouping specification

When multiple partitions exist, GTR parameters can be shared at three nested tiers:

- W groups (exchangeability): pool substitution counts $n_{ij}$ and time-in-state $\tau_j$ across partitions, update shared exchangeability matrix. Safe to share broadly (reflects chemistry, not selection).
- $\pi$ groups (equilibrium frequencies): pool row-summed counts and root state frequencies within groups. Share within functionally similar partitions (e.g., all 3rd codon positions).
- $\mu$ groups (overall rate): pool total counts, update shared rate. Typically per-partition (rates differ by function).

Nesting: W groups contain $\pi$ groups contain $\mu$ groups. A grouping specification defines which partitions share which parameters.

Default: share W broadly, $\pi$ within functionally similar groups, $\mu$ per-partition. Pseudocode for the generic inference function is in the auto-partitioning report.

**Sources:** `kb/reports/auto-partitioning.md` S7.2

### 4.2 Site-specific GTR as a single dense partition

Site-specific models (per-site $\pi^a$ and $\mu^a$, shared W) are one dense partition with a `SiteSpecific` model variant, not thousands of separate partitions. This preserves W sharing and traversal efficiency.

Proposed type: `enum GtrModel { Scalar(GTR), SiteSpecific(GTRSiteSpecific) }`.

The mathematical core is implemented and tested (4 golden master tests at 1e-10, 15 property tests). No production callers. Not wired into the partition system.

Design constraints: site-specific GTR is dense-only. Reject `--sparse` + `--site-specific-gtr`.

**Sources:** `kb/issues/N-gtr-site-specific-partition-integration.md`, `kb/reports/auto-partitioning.md` S7.3, R1, R2

### 4.3 Two-stage estimation (serialize and reload models)

Following IQ-TREE's PMSF workflow: estimate per-site preferences once on a guide tree, serialize to JSON (full model) or `.sitefreq` (IQ-TREE compatible), reload as fixed profiles for subsequent analyses. Decouples expensive estimation from repeated operations.

Two serialization proposals exist:

- **Scalar GTR:** JSON with `format: "treetime-gtr"`, fields for W (symmetric matrix), $\pi$ (frequencies), $\mu$ (rate scalar). Secondary exports: PAML triangular format, RAxML-NG model string. Enables `--model-file` and `--export-model` CLI flags.
- **Site-specific GTR:** JSON with `format: "treetime-site-specific-gtr"`, fields for W (shared matrix), $\pi$ (per-site: `n_states x seq_len`), $\mu$ (per-site: `seq_len`). Secondary export: IQ-TREE `.sitefreq`. Approximately 2 MB JSON for 30K-site nucleotide alignment.

**Sources:** `kb/proposals/model-serialization-scalar.md`, `kb/proposals/model-serialization-site-specific.md`, `kb/reports/auto-partitioning.md` S7.5

### 4.4 Edge-linked vs edge-unlinked branch lengths

Partitions sharing a tree can share branch lengths in two modes:

- Edge-proportional (scaled, RAxML-NG default): shared topology and branch lengths, per-partition rate scaler $\mu_p$. Branch length for partition $p$ on edge $e$: $b_{pe} = \mu_p \cdot b_e$.
- Edge-unlinked (independent): each partition has fully independent branch lengths. More parameters, less borrowing of signal.

v1 currently shares branch lengths at the graph level (implicitly edge-linked with $\mu$ from GTR). Design recommendation R9: adopt scaled partition as default.

**Sources:** `kb/proposals/config-file-multi-partition.md`, `kb/reports/auto-partitioning.md` S2.3, R9

## 5. Existing infrastructure

### 5.1 Partition type system

Three partition types with trait-based polymorphism:

| Type                      | Storage                                          | Use case                              |
| ------------------------- | ------------------------------------------------ | ------------------------------------- |
| `PartitionFitch`          | Per-node state sets                              | Parsimony reconstruction, compression |
| `PartitionMarginalDense`  | `Array2<f64>` (positions x states) per node/edge | Full probability matrices             |
| `PartitionMarginalSparse` | `BTreeMap<usize, VarPos>` per node/edge          | Variable positions only               |

Traits: `PartitionMarginalOps` (attach, backward, forward, extract, reconstruct), `PartitionCompressed` (sparse access), `HasLogLh`, `PartitionTimetreeOps` (contributions, topology reconciliation), `PartitionRerootOps`, composed into `PartitionTimetreeAll`.

Commands operate on `Vec<Arc<RwLock<dyn PartitionTimetreeAll>>>` (type alias `PartitionTimetreeAllVec`).

**Sources:** `kb/decisions/partition-system-architecture.md`

### 5.2 Existing branches

| Branch                                     | Content                                       | Status                 |
| ------------------------------------------ | --------------------------------------------- | ---------------------- |
| `worktree/feat/multi-segment-genome-input` | `--segment` flag, segment-aware FASTA loading | 9 commits, unmerged    |
| `worktree/feat/optimize-multi-alignment`   | Config parsing, multi-partition tests         | Unmerged               |
| `feat/abstract-partitions-enum`            | Partition abstraction via enum                | Historical, superseded |
| `feat/abstract-partitions-trait-objs`      | Partition abstraction via trait objects       | Became current design  |

## 6. Relationship to issue #2 (partition init duplication)

The duplication in partition init code across four commands is a symptom. The four commands do not create duplicated identical partitions. They each create independent partitions with different models, different dense/sparse decisions, and different post-init handling. The real problem is that partition creation is entangled with command orchestration instead of being a configured, composable layer.

The fix is not "extract shared init helper" but "implement partition configuration": commands receive configured partitions rather than constructing them inline. This requires:

1. A `PartitionSpec` type describing what partition to create (alignment, model, dense/sparse, index)
2. A `create_partitions(specs, graph) -> Vec<Partition>` function that builds partitions from specs
3. Configuration sources (CLI flags, config file, auto-detection) that produce `Vec<PartitionSpec>`
4. Commands that accept pre-built partitions

## 7. Cross-references

### Issues

- [N-optimize-multi-alignment-input](../issues/N-optimize-multi-alignment-input.md)
- [N-io-multi-segment-genome-input](../issues/N-io-multi-segment-genome-input.md)
- [N-representation-infer-dense-stub](../issues/N-representation-infer-dense-stub.md)
- [N-gtr-site-specific-partition-integration](../issues/N-gtr-site-specific-partition-integration.md)
- [M-core-partition-init-orchestration-duplication](../issues/M-core-partition-init-orchestration-duplication.md)

### Proposals

- [config-file-multi-partition](../proposals/config-file-multi-partition.md)
- [model-serialization-scalar](../proposals/model-serialization-scalar.md)
- [model-serialization-site-specific](../proposals/model-serialization-site-specific.md)

### Decisions

- [partition-system-architecture](../decisions/partition-system-architecture.md)

### Reports

- [auto-partitioning](auto-partitioning.md)

### Features

- [representation](../features/representation.md) - codon-position partitioning `[ ]`

### Design documents

- `_raw/optimize.md` - "alignment(s) (corresponding to partitions)"
- `_raw/sequence_evolution.md` - segment handling question
