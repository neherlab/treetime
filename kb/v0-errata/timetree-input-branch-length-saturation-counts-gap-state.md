# Input-mode branch-length variance counts the gap as a nucleotide state

## v0 location

`BranchLenInterpolator.__init__()` with `branch_length_mode == 'input'` computes the saturation term `p0 = 1 - sum(Pi^2)` from the GTR of the run [packages/legacy/treetime/treetime/branch_len_interpolator.py#L76](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L76). For long branches (input length of 0.05 substitutions per site and above), `p0` sets the variance of the Gaussian branch-length likelihood [packages/legacy/treetime/treetime/branch_len_interpolator.py#L87-L102](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L87-L102).

## Erratum

The default nucleotide alphabet `nuc` of v0 is `ACGT-`: the gap is a fifth state with its own equilibrium frequency. The default model without an alignment, `GTR.standard('jc', alphabet='nuc')` [packages/legacy/treetime/treetime/wrappers.py#L55-L56](../../packages/legacy/treetime/treetime/wrappers.py#L55-L56), has `Pi = 0.2` for each of the five states, so `p0 = 0.8`. With an inferred model, `Pi` includes the gap frequency of the alignment. The saturation term of nucleotide substitution is `1 - sum(Pi^2)` over the four bases, which is `0.75` for uniform frequencies.

## Evidence

- The v0 comment that derives the variance states that `p_0` is the distance of random sequences, "3/4 for nucleotides" [packages/legacy/treetime/treetime/branch_len_interpolator.py#L65-L75](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L65-L75); the code gives `0.8` with the default alphabet
- Under JC69 on four bases, two independent random sequences differ at a site with probability `1 - 4 * (1/4)^2 = 3/4`, the saturation level of the observed difference `p = p0 * (1 - exp(-l/p0))`. A gap is an indel, which no substitution process produces, so it is no state of this saturation
- v0 has the gap-free alphabet `nuc_nogap` (`ACGT`), whose JC69 model gives `p0 = 0.75`

## v0 impact

On long input branches, v0 overestimates `p0`, which underestimates the branch-length variance: by about 4% at 0.3 substitutions per site and about 25% at 1.5 substitutions per site with uniform frequencies and 1000 sites. Short branches (Poisson regime) do not use `p0`.

## v1 status

v1 computes `p0` over the states of its substitution model, which has four nucleotide states: JC69 without an alignment, the `--model` GTR with an alignment (`fn input_branch_model()` in [packages/treetime/src/timetree/inference/input_branch_length.rs](../../packages/treetime/src/timetree/inference/input_branch_length.rs)). The golden master of the input-mode likelihood captures v0 on the `nuc_nogap` alphabet, so the formula agrees with v0 and only the source of `Pi` differs.
