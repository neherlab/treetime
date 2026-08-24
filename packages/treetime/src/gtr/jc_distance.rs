#[cfg(test)]
mod __tests__;

/// Compute the Jukes-Cantor evolutionary distance from observed p-distance.
///
/// Under the Jukes-Cantor 1969 model (uniform equilibrium frequencies, equal
/// rates between all states), the observed proportion of differing sites $p$
/// (raw Hamming distance) between two sequences relates to the expected
/// evolutionary distance $d$ (substitutions per site) by inverting the
/// probability of observing a difference after time $t$ under a uniform
/// symmetric process on $k$ states:
///
/// $$ d = -\frac{k-1}{k} \, \ln\!\left(1 - \frac{k}{k-1}\, p\right) $$
///
/// Symbols:
///
/// - $p \in [0, (k-1)/k)$: observed p-distance (fraction of differing sites)
/// - $d \ge 0$: corrected evolutionary distance in expected substitutions per site
/// - $k \ge 2$: number of canonical states in the alphabet (4 for nucleotides,
///   20 for amino acids)
///
/// The raw p-distance always underestimates $d$ because repeated substitutions
/// at the same site may revert or mask earlier changes. For $k = 4$:
///
/// - $p = 0.01$: correction $\approx 0.7\%$
/// - $p = 0.10$: correction $\approx 7\%$
/// - $p = 0.25$: correction $\approx 22\%$
///
/// # Saturation
///
/// As $p$ approaches $p_{sat} = (k - 1) / k$, the corrected distance diverges
/// to infinity: at saturation, sequences are effectively random with respect
/// to each other under Jukes-Cantor and no finite $d$ can be inferred. To keep
/// the result finite and well-defined for downstream arithmetic (subtraction
/// from existing branch lengths, serialization), $p$ is clamped at
/// $p_{sat} \cdot (1 - \texttt{SATURATION\_MARGIN})$ before applying the
/// formula. The resulting cap is about 10 subs/site for $k = 4$ and 13 for
/// $k = 20$, large enough that any realistic child branch is consumed.
///
/// Negative $p$ is clamped to 0 defensively; callers should pass a valid ratio.
///
/// # Panics
///
/// Panics in debug builds if `n_states < 2`. The correction is undefined for
/// trivial alphabets; alphabets used throughout the crate always have at
/// least two canonical states.
///
/// # References
///
/// Jukes TH, Cantor CR (1969). Evolution of Protein Molecules. In: Munro HN
/// (ed.), Mammalian Protein Metabolism, vol. 3, pp. 21-132. Academic Press.
/// DOI: [10.1016/B978-1-4832-3211-9.50009-7](https://doi.org/10.1016/B978-1-4832-3211-9.50009-7)
pub fn jukes_cantor_distance(p: f64, n_states: usize) -> f64 {
  debug_assert!(
    n_states >= 2,
    "jukes_cantor_distance: n_states must be >= 2, got {n_states}"
  );

  let k = n_states as f64;
  let p_sat = (k - 1.0) / k;

  // Clamp into [0, p_sat * (1 - margin)) to avoid log(0), log(negative), or NaN
  let p_max = p_sat * (1.0 - SATURATION_MARGIN);
  let p_clamped = p.clamp(0.0, p_max);

  if p_clamped == 0.0 {
    return 0.0;
  }

  -(k - 1.0) / k * (1.0 - k / (k - 1.0) * p_clamped).ln()
}

/// Relative margin below the Jukes-Cantor saturation threshold.
///
/// Setting $p = p_{sat} \cdot (1 - \texttt{SATURATION\_MARGIN})$ bounds the
/// distance at roughly $-(k-1)/k \cdot \ln(\texttt{SATURATION\_MARGIN})$. The
/// value $10^{-6}$ gives a cap of about 10 substitutions per site for
/// nucleotide alphabets, which exceeds any realistic branch length while
/// avoiding the numerical issues of values closer to saturation.
const SATURATION_MARGIN: f64 = 1e-6;
