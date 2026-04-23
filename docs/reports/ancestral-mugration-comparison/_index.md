# Ancestral vs mugration: scientific and implementation comparison

A side-by-side comparison of the `ancestral` and `mugration` commands in both v0 (Python) and v1 (Rust). Covers the shared algorithmic core, divergent pipelines, data representations, GTR model handling, output formats, known issues, and architectural observations.

Companion to [Command relationships: ancestral, prune, optimize, timetree](../command-relationships/_index.md), which covers the four tree-refinement commands. That report does not include mugration. This report fills that gap.

## Core relationship

Both commands perform Felsenstein pruning (marginal maximum-likelihood ancestral reconstruction). Mugration is a specialization that applies the same algorithm to discrete traits (geographic location, host species) instead of biological sequences (nucleotides, amino acids).

The specialization changes three things: the alphabet, the sequence length, and the GTR model lifecycle.

| Aspect          | ancestral                                                       | mugration                                                 |
| :-------------- | :-------------------------------------------------------------- | :-------------------------------------------------------- |
| Input data      | FASTA alignment (hundreds to thousands of positions)            | CSV/TSV with one discrete trait per taxon                 |
| Alphabet        | Predefined nucleotide (4-5 states) or amino acid (20-22 states) | Custom: one state per unique trait value, arbitrary count |
| Sequence length | Full alignment length                                           | 1 (single pseudo-position encoding the trait)             |
| GTR lifecycle   | One-shot inference or named model                               | Iterative refinement (default 5 rounds)                   |

Everything else -- the tree structure, the message-passing direction, the Felsenstein recursion -- is the same mathematical operation.

## v0 (Python) pipeline comparison

### ancestral

Entry: `ancestral_reconstruction()` at [packages/legacy/treetime/treetime/wrappers.py#L588](../../../packages/legacy/treetime/treetime/wrappers.py#L588)

Pipeline:

1. Load tree and FASTA alignment into `TreeAnc(tree, aln=aln, gtr=gtr, ...)`
2. Call `treeanc.infer_ancestral_sequences('ml', infer_gtr=(gtr=='infer'), marginal=params.marginal, ...)`
   - `marginal` defaults to `False` (joint reconstruction)
   - `normalized_rate` defaults to `True` (mu forced to 1.0)
   - `pc` defaults to 5.0 (from `infer_gtr()` signature at [packages/legacy/treetime/treetime/treeanc.py#L1500](../../../packages/legacy/treetime/treetime/treeanc.py#L1500))
3. Write FASTA + Nexus + mutations TSV + Auspice JSON

Single call. No rate optimization. No iteration.

### mugration

Entry: `mugration()` at [packages/legacy/treetime/treetime/wrappers.py#L814](../../../packages/legacy/treetime/treetime/wrappers.py#L814), delegating to `reconstruct_discrete_traits()` at [packages/legacy/treetime/treetime/wrappers.py#L653](../../../packages/legacy/treetime/treetime/wrappers.py#L653)

Pipeline:

1. Map each unique trait value to a single uppercase letter (`"usa" -> "A"`, `"uk" -> "B"`, ...). Add a missing-data character as the next letter.
2. Build a custom GTR: `GTR.custom(pi=weights, W=ones((n,n)), alphabet=array(letters))`. Inject `profile_map[missing_char] = ones(n)`.
3. Create `TreeAnc(tree, gtr=mugration_GTR, ref='A', convert_upper=False, one_mutation=0.001, ...)`. Set `use_mutation_length = False`. Assign single-character pseudo-sequences: `{name: {0: letter}}`.
4. Initial reconstruction: `treeanc.infer_ancestral_sequences(method='ml', infer_gtr=True, marginal=True, normalized_rate=False, pc=pc, fixed_pi=weights, reconstruct_tip_states=True)`, then `treeanc.optimize_gtr_rate()`
5. Iterative refinement (5 rounds): `treeanc.infer_gtr(marginal=True, normalized_rate=False, pc=pc, fixed_pi=weights)` + `treeanc.optimize_gtr_rate()`
6. Optional sampling bias correction: `treeanc.gtr.mu *= correction`
7. Final reconstruction with refined, fixed GTR
8. Write GTR file + Nexus + optional confidence CSV

### v0 parameter differences

| Parameter                | ancestral                               | mugration                                               |
| :----------------------- | :-------------------------------------- | :------------------------------------------------------ |
| `marginal`               | User choice (default: `False` = joint)  | Always `True`                                           |
| `normalized_rate`        | `True` (mu = 1.0)                       | `False` (mu estimated)                                  |
| `pc`                     | 5.0 (default in `infer_gtr()`)          | 1.0 (default in `reconstruct_discrete_traits()`)        |
| `infer_gtr`              | Only if `--gtr infer`                   | Always `True`                                           |
| `reconstruct_tip_states` | User choice (default: `False`)          | Always `True`                                           |
| `optimize_gtr_rate()`    | Never called                            | Called after every GTR inference                        |
| `fixed_pi`               | `None` (unless VCF)                     | Weights from file or `None`                             |
| `one_mutation`           | `1/aln_length` (e.g. 0.001 for 1000 bp) | `1.0` (`full_length=1`, passed kwarg `0.001` is unused) |
| `store_compressed`       | Default `True`                          | `False`                                                 |
| `use_mutation_length`    | Default `False`                         | Explicitly set `False`                                  |

## v1 (Rust) pipeline comparison

### ancestral

Entry: `run_ancestral_reconstruction()` at [packages/treetime/src/commands/ancestral/run.rs#L35](../../../packages/treetime/src/commands/ancestral/run.rs#L35)

Pipeline (marginal dense path):

1. Read FASTA into `Vec<FastaRecord>`, read Newick into `GraphAncestral`
2. Create `PartitionMarginalDense` with dummy JC69 GTR
3. `initialize_marginal()` -- attach sequences, run backward + forward with dummy GTR
4. `get_gtr_dense()` -- infer real GTR from populated profiles
5. `update_marginal()` -- re-run backward + forward with inferred GTR
6. `ancestral_reconstruction_marginal()` -- extract MAP sequences
7. `annotate_branch_mutations()` -- populate mutation annotations
8. Write FASTA + Newick + Nexus + GTR JSON

Pipeline (marginal sparse path):

1. `compress_sequences()` -- Fitch backward + forward to build sparse representation
2. `get_gtr_sparse()` -- infer GTR from Fitch-compressed data
3. `update_marginal()` -- one marginal pass (no second pass)
4. Same reconstruction and annotation as dense

Pipeline (parsimony path):

1. `compress_sequences()` -- Fitch backward + forward
2. `ancestral_reconstruction_fitch()` -- extract sequences
3. No GTR needed

### mugration

Entry: `run_mugration()` at [packages/treetime/src/commands/mugration/run.rs#L258](../../../packages/treetime/src/commands/mugration/run.rs#L258), delegating to `parse_mugration_input()` + `execute_mugration()`

Pipeline:

1. Read Newick into `GraphAncestral`. Read traits from CSV/TSV into `BTreeMap<String, String>`. Optionally read weights.
2. Build `DiscreteStates` from unique trait values (sorted, deduplicated, missing marker excluded)
3. Compute pi: from weights (with mean-weight fallback for missing entries) or uniform. Apply pseudo-count smoothing if `--pc` is set.
4. Create `GTR::new(n_states, mu=1.0, W=None, pi=pi)` -- uniform rates
5. Create `PartitionDiscrete::new(0, gtr, discrete_states)`. Set `min_branch_length = 0.001`.
6. `attach_traits()` -- map leaf names to one-hot or uniform profiles, validate all leaves have traits and vice versa
7. `run_discrete_marginal()` -- backward + forward on `PartitionDiscrete`
8. `refine_gtr_iterative()` at [packages/treetime/src/commands/mugration/gtr_refinement.rs#L26](../../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L26):
   - Count transitions from marginal profiles via `count_transitions_discrete()`, then `infer_gtr_impl()`
   - `optimize_gtr_rate()` via custom Brent minimization
   - Repeat for `iterations` rounds (default 5)
   - Apply sampling bias correction
   - Final `run_discrete_marginal()` with refined model
9. Write Nexus + traits CSV + GTR JSON + optional confidence CSV

### v1 structural differences

| Aspect              | ancestral                                                                | mugration                                                                                       |
| :------------------ | :----------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------- |
| Alphabet type       | `Alphabet` (`#Alphabet`)                                                 | `DiscreteStates` (`#DiscreteStates`)                                                            |
| Partition types     | 3: `PartitionFitch`, `PartitionMarginalSparse`, `PartitionMarginalDense` | 1: `PartitionDiscrete`                                                                          |
| Profile shape       | 2D: `(positions, states)`                                                | 1D: `(states,)`                                                                                 |
| Trait dispatch      | `PartitionMarginalOps` trait                                             | Direct methods on `PartitionDiscrete`                                                           |
| Locking             | `Arc<RwLock<Partition>>`                                                 | `Arc<Mutex<&mut Partition>>`                                                                    |
| Dense/sparse toggle | Yes (`--dense` flag)                                                     | N/A (always dense over states)                                                                  |
| GTR init            | Named model or dummy JC69, then infer                                    | Uniform W + pi from data/weights                                                                |
| GTR refinement      | Single inference pass                                                    | Iterative (5 rounds) + Brent rate optimization                                                  |
| Rate optimization   | None                                                                     | `optimize_gtr_rate()` with Brent bracket check                                                  |
| Testability pattern | Monolithic `run_ancestral_reconstruction()`                              | Separated `parse_mugration_input()` + `execute_mugration()`                                     |
| Module count        | 4 (`args`, `fitch`, `marginal`, `run`)                                   | 7 (`args`, `comment_provider`, `discrete_marginal`, `gtr_refinement`, `input`, `output`, `run`) |
| Output              | FASTA + Newick + Nexus + GTR JSON                                        | Nexus + traits CSV + GTR JSON + optional confidence CSV                                         |

## Shared code

### Used by both commands

- `GraphAncestral` type alias ([packages/treetime/src/representation/payload/ancestral.rs#L15](../../../packages/treetime/src/representation/payload/ancestral.rs#L15)) -- both use the same graph with `NodeAncestral` / `EdgeAncestral` payloads
- `nwk_read_file()` -- tree input
- `GTR` struct and `GTR::new()` constructor ([packages/treetime/src/gtr/gtr.rs](../../../packages/treetime/src/gtr/gtr.rs))
- `infer_gtr_impl()` from `gtr/infer_gtr/common` -- core GTR inference algorithm, called by ancestral (via `infer_gtr_dense` / `infer_gtr_sparse`) and by mugration (via `gtr_refinement.rs`)
- `get_branch_mutation_matrix()` and `accumulate_mutation_counts()` from `gtr/infer_gtr/dense` -- reused by mugration's `count_transitions_discrete()` ([packages/treetime/src/commands/mugration/gtr_refinement.rs#L128](../../../packages/treetime/src/commands/mugration/gtr_refinement.rs#L128))
- `softmax_with_log_norm()` from `partition/marginal_helpers` -- used by `PartitionDiscrete::process_node_backward()` ([packages/treetime/src/representation/partition/discrete.rs#L82](../../../packages/treetime/src/representation/partition/discrete.rs#L82))

### Ancestral only

- `PartitionFitch`, `PartitionMarginalDense`, `PartitionMarginalSparse` partition types
- Fitch parsimony algorithm ([packages/treetime/src/commands/ancestral/fitch.rs](../../../packages/treetime/src/commands/ancestral/fitch.rs))
- Generic marginal backward/forward via `PartitionMarginalOps` trait ([packages/treetime/src/commands/ancestral/marginal.rs](../../../packages/treetime/src/commands/ancestral/marginal.rs))
- `annotate_branch_mutations()` for mutation annotation on graph nodes
- `Alphabet`, `Seq`, `StateSet`, `BitSet128` sequence primitives

### Mugration only

- `PartitionDiscrete` with self-contained `process_node_backward()` / `process_node_forward()` ([packages/treetime/src/representation/partition/discrete.rs](../../../packages/treetime/src/representation/partition/discrete.rs))
- `DiscreteStates` instead of `Alphabet` ([packages/treetime/src/representation/discrete_states.rs](../../../packages/treetime/src/representation/discrete_states.rs))
- `DiscreteNodeData` / `DiscreteEdgeData` payload types ([packages/treetime/src/representation/payload/discrete.rs](../../../packages/treetime/src/representation/payload/discrete.rs))
- Iterative GTR refinement with custom Brent optimizer ([packages/treetime/src/commands/mugration/gtr_refinement.rs](../../../packages/treetime/src/commands/mugration/gtr_refinement.rs))
- Sampling bias correction
- Weight-based pi computation with pseudo-count smoothing
- `PartitionCommentProvider` for Nexus trait annotation ([packages/treetime/src/commands/mugration/comment_provider.rs](../../../packages/treetime/src/commands/mugration/comment_provider.rs))
- Structured I/O types: `MugrationInput`, `MugrationResult`, `MugrationGtrOutput`, `MugrationTraitsOutput`, `MugrationConfidenceOutput`

## GTR model lifecycle

The GTR model lifecycle is the single largest algorithmic difference between the two commands. This is consistent across v0 and v1.

### ancestral GTR lifecycle

```
named model (JC69, K80, ...) or dummy JC69
    |
    v
(optional) single infer_gtr pass from Fitch/marginal profiles
    |
    v
fixed GTR for marginal reconstruction
```

- `normalized_rate=True` in v0: mu forced to 1.0 after inference
- v1: no rate optimization, GTR used as-is after inference

### mugration GTR lifecycle

```
uniform W + pi (from weights or uniform)
    |
    v
initial marginal reconstruction
    |
    v
infer GTR from marginal profiles + optimize mu (Brent)
    |
    +---> repeat 5 times: re-infer GTR + re-optimize mu
    |
    v
(optional) multiply mu by sampling bias correction
    |
    v
final marginal reconstruction with refined GTR
```

- `normalized_rate=False` in v0: mu estimated freely
- v1: custom Brent minimization of `-log_lh` over `sqrt(mu)`, with bracket validation

### Why mugration needs iteration

Mugration operates on a single "position" with $N$ discrete states. Small shifts in $\pi$ (equilibrium frequencies) or $W$ (rate matrix) substantially change posteriors at ambiguous internal nodes, because there is no averaging across hundreds of alignment columns. Ancestral reconstruction aggregates over many sites, so expected substitution counts constrain the GTR model more tightly from the first inference pass.

A proposal to add iterative GTR refinement to ancestral exists ([docs/port-proposals/ancestral-iterative-gtr-refinement.md](../../../docs/port-proposals/ancestral-iterative-gtr-refinement.md)) but is not accepted. That proposal notes iterative refinement would help ancestral in specific cases: short alignments, strong compositional bias, large alphabets, weakly resolved branches, poor initialization model.

## Message passing implementation

Both commands implement the same Felsenstein pruning algorithm (backward: leaves to root, forward: root to leaves). The implementations diverge in profile shape and trait dispatch.

### ancestral (via `PartitionMarginalOps` trait)

- Profiles are 2D arrays: `(n_positions, n_states)`
- Backward/forward dispatch through `PartitionMarginalOps` trait, implemented by both `PartitionMarginalDense` and `PartitionMarginalSparse`
- Generic functions `marginal_backward()` / `marginal_forward()` in [packages/treetime/src/commands/ancestral/marginal.rs](../../../packages/treetime/src/commands/ancestral/marginal.rs) iterate over partitions behind `Arc<RwLock<P>>`
- Transition matrix $e^{Qt}$ applied per-position

### mugration (direct on `PartitionDiscrete`)

- Profiles are 1D arrays: `(n_states,)`
- Backward/forward implemented directly as methods on `PartitionDiscrete` ([packages/treetime/src/representation/partition/discrete.rs#L59-L213](../../../packages/treetime/src/representation/partition/discrete.rs#L59-L213))
- Does not implement `PartitionMarginalOps`
- Traversal uses `Arc<Mutex<&mut Partition>>` instead of `Arc<RwLock<Partition>>`
- Single transition matrix $e^{Qt}$ per edge (one "position")

### Backward pass differences

| Step           | ancestral (dense)                                                                            | mugration                                                       |
| :------------- | :------------------------------------------------------------------------------------------- | :-------------------------------------------------------------- |
| Leaf profile   | 2D one-hot from sequence per position                                                        | 1D one-hot from observed trait (or uniform if missing)          |
| Internal node  | Product of child messages per position, normalize                                            | Product of child messages in log space, `softmax_with_log_norm` |
| Root weighting | Per-position: `profile * pi`                                                                 | Single: `profile * pi`                                          |
| Edge message   | $\text{msg\_from\_child}[j] = \sum_i \text{msg\_to\_parent}[i] \cdot P_{ij}(t)$ per position | Same formula, 1D                                                |

### Forward pass differences

| Step                    | ancestral (dense)                                             | mugration                                                      |
| :---------------------- | :------------------------------------------------------------ | :------------------------------------------------------------- |
| Parent-to-child message | `msg_to_child = parent_profile / msg_from_child` per position | Same formula, 1D, with `f64::MIN_POSITIVE` clamping on divisor |
| Propagation             | Through $e^{Qt}$ per position                                 | Through $e^{Qt}$, single                                       |
| Combination             | Product of propagated parent message and subtree message      | Same                                                           |

## Data representations

### Node data

| Field          | ancestral dense (`DenseNodePartition`)     | mugration (`DiscreteNodeData`)          |
| :------------- | :----------------------------------------- | :-------------------------------------- |
| Observed       | Original sequence (via `Seq`)              | `observed: Option<usize>` (state index) |
| Profile        | `Array2<f64>` -- `(n_positions, n_states)` | `Array1<f64>` -- `(n_states,)`          |
| Log-likelihood | Accumulated in node                        | `log_lh: f64`                           |

### Edge data

| Field            | ancestral dense (`DenseEdgePartition`) | mugration (`DiscreteEdgeData`) |
| :--------------- | :------------------------------------- | :----------------------------- |
| `msg_to_parent`  | `Array2<f64>`                          | `Array1<f64>`                  |
| `msg_from_child` | `Array2<f64>`                          | `Array1<f64>`                  |
| `msg_to_child`   | `Array2<f64>`                          | `Array1<f64>`                  |
| Log-likelihood   | Per-edge                               | `log_lh: f64`                  |

### Missing data handling

- **ancestral**: Ambiguous characters (N, gaps) produce uniform-like profiles at affected positions via `Alphabet.profile_map`. Terminal gap handling controlled by `--keep-overhangs`.
- **mugration**: Missing data marker (default `"?"`) produces `DiscreteNodeData::missing(n_states)` = uniform profile `1/n_states` across all states. `DiscreteStates` excludes the missing marker from the state set entirely.

## Stale-messages bug in iterative refinement

Both v0 and v1 share a subtle correctness issue in the iterative GTR refinement loop. The loop re-estimates the GTR model from marginal profiles but does **not** re-run a full forward-backward pass between iterations. This means:

- Forward messages from the initial reconstruction become stale after the GTR changes
- `count_transitions_discrete()` / v0's `infer_gtr()` use these stale forward messages to compute expected transition counts
- The resulting EM does not satisfy the monotonicity guarantee of a proper EM algorithm

A proposal to fix this exists: [docs/port-proposals/mugration-full-reconstruction-per-iteration.md](../../../docs/port-proposals/mugration-full-reconstruction-per-iteration.md). The fix would insert a full `run_discrete_marginal()` call before each `count_transitions_discrete()` in the iterative loop. Expected impact: correct EM iterations, different fixed point from both current v1 and v0, ~2x computational cost per iteration.

## Known issues and intentional changes

### Ancestral known issues (7 open)

| ID                                         | Issue                                                                                       | Impact                                               |
| :----------------------------------------- | :------------------------------------------------------------------------------------------ | :--------------------------------------------------- |
| `H-ancestral-joint-default-panics`         | `--method-anc` defaults to `Joint` which panics (`unimplemented!()`)                        | Default invocation crashes                           |
| `M-ancestral-dense-sparse-divergence`      | ~2.5% of random gap-free GTR configs produce different likelihoods between dense and sparse | Affects sparse path only                             |
| `M-ancestral-marginal-probability-space`   | Forward pass uses plain probability space (v1) vs neg-log space (v0)                        | Floating-point divergence, GM tolerance 1e-6 to 1e-7 |
| `M-ancestral-sparse-root-invariance`       | Sparse marginal violates Felsenstein pulley principle by ~0.09                              | Affects sparse path only                             |
| `M-ancestral-sparse-alphabet-mismatch`     | Variable-site classification uses `is_ambiguous()` vs v0's `alphabet_gapN`                  | Contributes to dense-sparse divergence               |
| `M-ancestral-stdin-fasta-truncation`       | Stdin FASTA reads one record, silently truncates multi-record alignments                    | Data loss on piped input                             |
| `N-ancestral-sparse-remove-insert-pattern` | Remove/insert pattern in sparse passes                                                      | Code smell                                           |
| ~~`N-dense-normalize-inplace-zero-row`~~   | ~~NaN for all-zero probability rows~~ (fixed: returns uniform distribution)                 | ~~Edge case guard~~                                  |

None of these affect mugration (mugration uses `PartitionDiscrete`, not the sequence-based partitions).

### Mugration known issues (1 open)

| ID                          | Issue                                                                                         | Impact                             |
| :-------------------------- | :-------------------------------------------------------------------------------------------- | :--------------------------------- |
| `M-mugration-iterative-gtr` | 5/7 datasets diverge from v0 due to D1 (pseudo-count smoothing) and D2 (root state filtering) | Intentional improvements, not bugs |

### Ancestral intentional changes (2)

| Change                                     | Description                                                                                         |
| :----------------------------------------- | :-------------------------------------------------------------------------------------------------- |
| `ancestral-joint-reconstruction-removed`   | Joint ML removed -- statistically inconsistent per Mossel et al. 2009. `--method-anc joint` panics. |
| `ancestral-fitch-deterministic-root-state` | v1 uses `BitSet128::get_one()` (lowest ASCII). v0 uses random selection.                            |

### Mugration intentional changes (2)

| Change                              | Description                                                                                                                            |
| :---------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------- |
| `mugration-pseudo-count-initial-pi` | Pseudo-count smoothing on initial pi. Smoother Bayesian prior (Dirichlet).                                                             |
| `mugration-root-state-filtering`    | Root state excluded from pi estimation when posterior is near-uniform (`max_prob <= 1/n_states + 1e-10`). Eliminates state-order bias. |

### Shared intentional change (1)

| Change                                   | Description                                                                                                                       |
| :--------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------- |
| `gtr-uninformative-root-state-filtering` | Dense GTR inference skips uninformative positions when computing `root_state`. Mugration's root filtering is the discrete analog. |

### Feature parity with v0

- **Ancestral**: 10/16 features implemented, 6 missing. The missing features are CLI flags that are parsed but not wired: `--keep-overhangs`, `--zero-based`, `--report-ambiguous`, `--seed`, `--gtr-params`, `--aa`, `--vcf-reference`, `--aln`.
- **Mugration**: 10/10 features implemented, 0 missing.

## Architectural observations

### Mugration does not participate in the command-relationships hierarchy

The [command-relationships report](../command-relationships/_index.md) identifies `ancestral` as the shared foundation for `prune`, `optimize`, and `timetree`. Mugration is outside this hierarchy entirely:

- It does not share partition types with any other command
- It does not share the `PartitionMarginalOps` trait
- It does not feed into optimize or timetree loops
- It does not use Fitch compression

Mugration is a **sibling** of the ancestral-prune-optimize-timetree family, not a member. It shares the GTR inference core (`infer_gtr_impl`, `get_branch_mutation_matrix`, `accumulate_mutation_counts`) and the graph structure (`GraphAncestral`), but its message-passing, data representation, and iterative refinement are self-contained.

### `PartitionDiscrete` reimplements rather than specializes

`PartitionDiscrete` implements its own `process_node_backward()` / `process_node_forward()` methods directly, without implementing the `PartitionMarginalOps` trait. This is a pragmatic choice -- 1D `Array1<f64>` profiles are structurally different from 2D `Array2<f64>` profiles, and forcing them through the same trait would require either generics over dimensionality or reshape overhead.

The cost is code duplication: the Felsenstein recursion (product of child messages, normalize, propagate through transition matrix, divide parent profile by child message) is written twice with minor differences. Changes to one implementation must be manually synchronized to the other.

### Mugration's separated input/output is a better pattern

Mugration separates `parse_mugration_input()` (file I/O to `MugrationInput`) from `execute_mugration()` (`MugrationInput` to `MugrationResult`). This makes the execution core testable without file I/O.

Ancestral has a monolithic `run_ancestral_reconstruction()` that interleaves file reading, partition creation, reconstruction, and output writing. This pattern is shared by other tree-refinement commands.

### Locking strategy differs

Ancestral uses `Arc<RwLock<Partition>>` -- multiple partitions can be read concurrently while one is being written. Mugration uses `Arc<Mutex<&mut Partition>>` -- exclusive access only. Since mugration has exactly one partition and the traversal is inherently sequential (each node depends on its parent/children), the simpler `Mutex` is sufficient.
