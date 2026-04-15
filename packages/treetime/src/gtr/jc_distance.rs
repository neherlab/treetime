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

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use rstest::rstest;

  #[test]
  fn test_jukes_cantor_distance_zero_p_returns_zero() {
    // Bit-exact: the zero-input branch returns the 0.0 literal
    assert_abs_diff_eq!(0.0, jukes_cantor_distance(0.0, 4), epsilon = 0.0);
    assert_abs_diff_eq!(0.0, jukes_cantor_distance(0.0, 20), epsilon = 0.0);
  }

  #[test]
  fn test_jukes_cantor_distance_negative_p_clamped_to_zero() {
    // Defensive: negative inputs clamp to zero rather than producing NaN
    assert_abs_diff_eq!(0.0, jukes_cantor_distance(-0.1, 4), epsilon = 0.0);
    assert_abs_diff_eq!(0.0, jukes_cantor_distance(-1.0, 4), epsilon = 0.0);
  }

  #[rustfmt::skip]
  #[rstest]
  // Reference values computed in Python with the bit-exact formula
  //   d = -(k-1)/k * log(1 - k/(k-1) * p)
  // using double-precision math.log(). These match Rust f64 semantics.
  #[case::nuc_small_p(      (0.01,  4), 0.010067265249105496)]
  #[case::nuc_issue_example((0.10,  4), 0.10732563273050497 )]
  #[case::nuc_large_p(      (0.25,  4), 0.30409883108112323 )]
  #[case::nuc_half_sat(     (0.375, 4), 0.5198603854199589  )]
  #[case::aa_small_p(       (0.05, 20), 0.05136386020676192 )]
  #[case::aa_medium_p(      (0.20, 20), 0.22456933916101887 )]
  #[trace]
  fn test_jukes_cantor_distance_known_values(
    #[case] (p, n_states): (f64, usize),
    #[case] expected: f64,
  ) {
    let actual = jukes_cantor_distance(p, n_states);
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-15);
  }

  #[rustfmt::skip]
  #[rstest]
  // Saturation threshold (k-1)/k: finite cap, never NaN or infinity
  #[case::nuc_at_saturation(    (0.75,  4))]
  #[case::nuc_beyond_saturation((0.80,  4))]
  #[case::nuc_p_one(            (1.00,  4))]
  #[case::aa_at_saturation(     (0.95, 20))]
  #[case::aa_beyond_saturation( (0.99, 20))]
  #[trace]
  fn test_jukes_cantor_distance_saturation_is_finite(#[case] (p, n_states): (f64, usize)) {
    let d = jukes_cantor_distance(p, n_states);
    assert!(d.is_finite(), "expected finite distance at saturation, got {d}");
    assert!(d > 0.0, "expected positive distance at saturation, got {d}");
  }

  #[test]
  fn test_jukes_cantor_distance_saturation_cap_order_of_magnitude() {
    // For k=4, SATURATION_MARGIN=1e-6 gives d_max ~ 0.75 * ln(1e6) ~ 10.36
    let d_max_nuc = jukes_cantor_distance(1.0, 4);
    assert!((10.0..11.0).contains(&d_max_nuc), "d_max nuc = {d_max_nuc}");
    // For k=20, d_max ~ 0.95 * ln(1e6) ~ 13.13
    let d_max_aa = jukes_cantor_distance(1.0, 20);
    assert!((12.0..14.0).contains(&d_max_aa), "d_max aa = {d_max_aa}");
  }

  #[test]
  fn test_jukes_cantor_distance_always_at_least_p() {
    // JC69 correction never underestimates: d >= p for all valid p
    for n in [2_usize, 4, 20] {
      let k = n as f64;
      let p_sat = (k - 1.0) / k;
      for i in 0..100 {
        let p = p_sat * 0.999 * f64::from(i) / 100.0;
        let d = jukes_cantor_distance(p, n);
        assert!(d >= p - 1e-15, "d={d} < p={p} for n={n}");
      }
    }
  }

  #[test]
  fn test_jukes_cantor_distance_monotonic_in_p() {
    // d(p) is strictly increasing on [0, p_sat)
    for n in [2_usize, 4, 20] {
      let k = n as f64;
      let p_sat = (k - 1.0) / k;
      let mut prev = 0.0;
      for i in 0..1000 {
        let p = p_sat * 0.999 * f64::from(i) / 1000.0;
        let d = jukes_cantor_distance(p, n);
        assert!(d >= prev, "non-monotonic at p={p} n={n}: {d} < {prev}");
        prev = d;
      }
    }
  }

  #[test]
  fn test_jukes_cantor_distance_small_p_approaches_p() {
    // For p -> 0, d(p) -> p (Taylor expansion: d ~ p + p^2/(2(k-1)/k) + ...)
    for p in [1e-6, 1e-5, 1e-4] {
      let d = jukes_cantor_distance(p, 4);
      let rel_err = (d - p).abs() / p;
      assert!(rel_err < 1e-3, "d={d} not close to p={p}: rel_err={rel_err}");
    }
  }

  #[test]
  fn test_jukes_cantor_distance_issue_documented_error() {
    // Issue doc: raw p underestimates d by 7% at p=0.1, 18% at p=0.25
    // Our exact corrections: ~7.3% at p=0.1, ~21.6% at p=0.25
    let d_0_1 = jukes_cantor_distance(0.1, 4);
    let rel_correction_0_1 = (d_0_1 - 0.1) / 0.1;
    assert!(
      (0.06..0.08).contains(&rel_correction_0_1),
      "expected ~7% correction at p=0.1, got {rel_correction_0_1}"
    );

    let d_0_25 = jukes_cantor_distance(0.25, 4);
    let rel_correction_0_25 = (d_0_25 - 0.25) / 0.25;
    assert!(
      (0.18..0.24).contains(&rel_correction_0_25),
      "expected ~18-22% correction at p=0.25, got {rel_correction_0_25}"
    );
  }
}
