# Dense initial GTR inference uses Fitch parsimony counts

## Change

Dense `--model=infer` now infers GTR parameters from Fitch parsimony substitution counts (`infer_gtr_fitch`) instead of from marginal posterior fractional counts (`infer_gtr_dense`). Both dense and sparse partitions use the same Fitch-based inference path for initial GTR estimation.

## v0 behavior

v0's `TreeAnc.infer_gtr(marginal=True)` infers GTR from branch joint distributions computed during marginal ancestral reconstruction. This requires a full marginal pass before GTR inference can run, creating a chicken-and-egg problem: partitions need a GTR to compute transition probabilities, but GTR inference needs populated partitions.

v0 works around this by initializing with a JC69 placeholder, running marginal reconstruction, inferring GTR from the (JC69-biased) profiles, then rerunning marginal with the inferred GTR. The second marginal pass corrects the JC69 bias but adds computational cost.

## v1 behavior

v1 runs Fitch parsimony compression (which does not require a GTR model) to produce integer substitution counts, then infers GTR from those counts. The inferred GTR is available before any marginal pass, eliminating the placeholder and double-pass workaround.

The `PartitionFitch` typestate lifecycle enforces this ordering: `compress` -> `infer_gtr` -> `into_marginal_dense(gtr)`.

## Measured deviations

Dense marginal-posterior GTR and Fitch parsimony GTR produce different parameters on real datasets (measured via `test_contract_dense_sparse_real.rs`):

| Dataset          | pi cosine   | W rel Frobenius | mu rel diff |
| ---------------- | ----------- | --------------- | ----------- |
| flu/h3n2/20      | 0.999998594 | 0.0173          | 0.0040      |
| ebola/20         | 0.999999994 | 0.0180          | 0.0236      |
| rsv/a/20         | 0.999991113 | 0.0070          | 0.0229      |
| dengue/20        | 0.999786491 | 0.0080          | 0.0413      |
| tb/20            | 0.999970946 | 0.0109          | 0.0247      |
| lassa/L/50       | 0.998686724 | 0.1015          | 0.0196      |
| mpox/clade-ii/20 | 0.999999989 | 0.0858          | 0.1126      |

Equilibrium frequencies (pi) are nearly identical (cosine > 0.997). Exchangeability matrix (W) differs up to ~10% relative Frobenius on datasets with high ambiguity (lassa, mpox). Rate (mu) differs up to ~11% on mpox (long genome, few mutations).

## Rationale

- Fitch counts are GTR-independent, eliminating the circular dependency
- Dense and sparse initial GTR are now identical, removing a source of discrepancy between the two modes
- Single marginal pass instead of double reduces computation
- The dense marginal-posterior path (`infer_gtr_dense`) is retained for future iterative GTR refinement

## Related

- [Dense partitions lack Fitch compression](../port-known-issues/H-dense-with-fitch-compression.md) - broader context for Fitch integration into dense mode
- [Sparse fixed-position scalar rate approximation](sparse-fixed-position-scalar-rate-approximation.md) - related sparse-specific deviation
