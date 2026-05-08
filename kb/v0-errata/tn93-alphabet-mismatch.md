# TN93 model ignores kappa parameters due to alphabet mismatch

## v0 location

`TN93()` (`#TN93`) [packages/legacy/treetime/treetime/nuc_models.py#L193-L240](../../packages/legacy/treetime/treetime/nuc_models.py#L193-L240)

## Erratum

`TN93()` constructs a 4x4 rate matrix W with kappa1/kappa2 transition/transversion parameters, validates pi against the 4-state `nuc_nogap` alphabet, but passes the 5-state `nuc` alphabet (A,C,G,T,-) to the GTR constructor. `GTR.assign_rates()` (`#assign_rates`) detects the dimension mismatch (4x4 matrix vs 5-state alphabet) and silently replaces W with a uniform 5x5 rate matrix. TN93 behaves identically to JC69.

## Evidence

- [packages/legacy/treetime/treetime/nuc_models.py#L224-L227](../../packages/legacy/treetime/treetime/nuc_models.py#L224-L227): creates 4x4 W matrix with kappa1/kappa2
- [packages/legacy/treetime/treetime/nuc_models.py#L230](../../packages/legacy/treetime/treetime/nuc_models.py#L230): validates pi against `alphabets['nuc_nogap']` (4 states) - correct
- [packages/legacy/treetime/treetime/nuc_models.py#L238](../../packages/legacy/treetime/treetime/nuc_models.py#L238): passes `alphabet=alphabets['nuc']` (5 states) to GTR constructor - wrong
- `GTR.assign_rates()` (`#assign_rates`) at [packages/legacy/treetime/treetime/gtr.py#L266-L271](../../packages/legacy/treetime/treetime/gtr.py#L266-L271): logs a warning at verbosity=4 (rarely seen) and replaces W with uniform rates
- All other nucleotide models (JC69, K2P, F81, HKY85) use `nuc_nogap` - TN93 is the outlier

## v0 impact

TN93 is non-functional. Any analysis specifying `--gtr=tn93` runs with JC69 rates.

## v1 status

v1's GTR system uses explicit alphabet types with compile-time dimension checks. The mismatch category cannot occur.

## References

- Tamura & Nei (1993). "Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees." Mol Biol Evol, 10(3):512-526. doi:10.1093/oxfordjournals.molbev.a040023
