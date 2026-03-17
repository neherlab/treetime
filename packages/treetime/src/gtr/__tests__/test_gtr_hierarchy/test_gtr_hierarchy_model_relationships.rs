#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{F81Params, HKY85Params, JC69Params, K80Params, TN93Params, f81, hky85, jc69, k80, tn93};
  use crate::gtr::gtr::{GTR, GTRParams};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::array;

  // Model Hierarchy Tests
  //
  // The nucleotide substitution models form a nested hierarchy where simpler
  // models are special cases of more complex ones:
  //
  //     JC69 -----> K80 -----> HKY85 -----> TN93 -----> GTR
  //       |                      ^
  //       +-------> F81 ---------+
  //
  // These tests verify that setting appropriate parameters in a more complex
  // model yields identical Q matrices to simpler models.
  //
  // NOTE: This implementation uses a transposed Q convention where Q[i,j] is
  // the rate from state j to state i. The tests compare normalized Q matrices
  // after mu normalization.

  /// JC69 is K80 with equal transition/transversion rates (kappa = 1).
  /// K80 distinguishes transitions (A<->G, C<->T) from transversions.
  /// When kappa = 1, this distinction vanishes, recovering JC69.
  /// Both have uniform equilibrium frequencies.
  #[test]
  fn test_gtr_jc69_equals_k80_kappa_1() -> Result<(), Report> {
    let jc = jc69(JC69Params {
      mu: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let k = k80(K80Params {
      mu: 1.0,
      kappa: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    // Compare normalized Q matrices
    let q_jc = jc.Q();
    let q_k = k.Q();
    assert_abs_diff_eq!(q_jc, q_k, epsilon = 1e-14);

    // Also verify pi is the same
    assert_abs_diff_eq!(jc.pi, k.pi, epsilon = 1e-14);

    Ok(())
  }

  /// JC69 is F81 with uniform equilibrium frequencies (pi = [0.25, 0.25, 0.25, 0.25]).
  /// F81 allows non-uniform pi but keeps all exchangeabilities equal.
  /// With uniform pi, the Q matrix becomes symmetric, recovering JC69.
  #[test]
  fn test_gtr_jc69_equals_f81_uniform_pi() -> Result<(), Report> {
    let jc = jc69(JC69Params {
      mu: 1.0,
      alphabet: AlphabetName::Nuc,
    })?;

    let f = f81(F81Params {
      mu: 1.0,
      pi: Some(array![0.25, 0.25, 0.25, 0.25]),
      alphabet: AlphabetName::Nuc,
    })?;

    // Compare normalized Q matrices
    let q_jc = jc.Q();
    let q_f = f.Q();
    assert_abs_diff_eq!(q_jc, q_f, epsilon = 1e-14);

    // Also verify pi is the same
    assert_abs_diff_eq!(jc.pi, f.pi, epsilon = 1e-14);

    Ok(())
  }

  /// K80 is HKY85 with uniform equilibrium frequencies.
  /// HKY85 combines the transition/transversion distinction of K80 with
  /// the non-uniform frequencies of F81. Restricting to uniform pi
  /// removes the frequency effect, recovering K80.
  #[test]
  fn test_gtr_k80_equals_hky85_uniform_pi() -> Result<(), Report> {
    let kappa = 2.5;

    let k = k80(K80Params {
      mu: 1.0,
      kappa,
      alphabet: AlphabetName::Nuc,
    })?;

    let h = hky85(HKY85Params {
      mu: 1.0,
      kappa,
      pi: Some(array![0.25, 0.25, 0.25, 0.25]),
      alphabet: AlphabetName::Nuc,
    })?;

    // Compare normalized Q matrices
    let q_k = k.Q();
    let q_h = h.Q();
    assert_abs_diff_eq!(q_k, q_h, epsilon = 1e-14);

    // Also verify pi is the same
    assert_abs_diff_eq!(k.pi, h.pi, epsilon = 1e-14);

    Ok(())
  }

  /// F81 is HKY85 with equal transition/transversion rates (kappa = 1).
  /// HKY85 with kappa = 1 treats all substitutions equally (like F81),
  /// while still allowing non-uniform equilibrium frequencies.
  #[test]
  fn test_gtr_f81_equals_hky85_kappa_1() -> Result<(), Report> {
    let pi = array![0.1, 0.2, 0.3, 0.4];

    let f = f81(F81Params {
      mu: 1.0,
      pi: Some(pi.clone()),
      alphabet: AlphabetName::Nuc,
    })?;

    let h = hky85(HKY85Params {
      mu: 1.0,
      kappa: 1.0,
      pi: Some(pi),
      alphabet: AlphabetName::Nuc,
    })?;

    // Compare normalized Q matrices
    let q_f = f.Q();
    let q_h = h.Q();
    assert_abs_diff_eq!(q_f, q_h, epsilon = 1e-14);

    // Also verify pi is the same
    assert_abs_diff_eq!(f.pi, h.pi, epsilon = 1e-14);

    Ok(())
  }

  /// HKY85 is TN93 with equal purine/pyrimidine transition rates (kappa1 = kappa2).
  /// TN93 distinguishes A<->G transitions from C<->T transitions.
  /// When both have the same rate, this reduces to HKY85's single kappa.
  ///
  /// NOTE: TN93 uses kappa1 for transversions and kappa2 for C<->T transitions,
  /// with A<->G as reference (rate=1). So TN93(kappa1=k, kappa2=1) gives
  /// the same relative rates as HKY85(kappa=1/k) when k is the tv/ts ratio.
  ///
  /// To match HKY85 exactly, we set TN93 with:
  /// - kappa1 = 1 (transversions same as A<->G reference)
  /// - kappa2 = 1 (C<->T same as A<->G)
  ///
  /// This gives all rates equal, like F81 with uniform pi = JC69.
  ///
  /// For non-trivial comparison, we use:
  /// - HKY85 with kappa (ts/tv ratio)
  /// - TN93 with kappa1=1/kappa (tv rate), kappa2=1 (C<->T same as A<->G)
  #[test]
  fn test_gtr_hky85_equals_tn93_equal_transitions() -> Result<(), Report> {
    let pi = array![0.1, 0.2, 0.3, 0.4];

    // HKY85: kappa is the transition/transversion ratio
    // In HKY85, W matrix has transitions (A<->G, C<->T) = kappa, transversions = 1
    //
    // TN93: kappa1 is transversion rate, kappa2 is C<->T rate, A<->G = 1 (reference)
    // For equal transitions: kappa2 = 1 (C<->T same as A<->G)
    // To get ts/tv ratio of kappa: need transversions = 1/kappa relative to transitions = 1
    // But TN93 has transversions = kappa1 and A<->G = 1
    // So for HKY85(kappa) ~ TN93(kappa1, kappa2), we need:
    //   - Both transitions equal: kappa2 = 1 (C<->T = A<->G = 1)
    //   - ts/tv ratio = 1/kappa1 = kappa => kappa1 = 1/kappa... but that's inverted
    //
    // Actually, looking at the W matrices:
    // HKY85 W: off-diag = 1 except transitions = kappa
    // TN93 W: off-diag = kappa1 except A<->G = 1, C<->T = kappa2
    //
    // For equivalence with kappa2=1 (equal transitions):
    // HKY85(kappa) matches TN93(1, 1) only when kappa=1 (JC69-like)
    //
    // For general equivalence where HKY85 kappa maps to TN93:
    // We need the W matrices to be proportional (same after normalization).
    // HKY85 W: transversions=1, transitions=kappa
    // TN93 W with kappa2=1: A<->G=1, C<->T=1, transversions=kappa1
    //
    // These match when kappa1=1 and kappa=1, or when HKY85 has kappa=k and
    // TN93 has kappa1=1/k, kappa2=1 (but W is normalized, so scale doesn't matter).
    //
    // Let's set kappa=2 for HKY85, and for TN93 use kappa1=0.5, kappa2=1.
    // After normalization, these should match.
    let h = hky85(HKY85Params {
      mu: 1.0,
      kappa: 2.0,
      pi: Some(pi.clone()),
      alphabet: AlphabetName::Nuc,
    })?;

    // TN93 with equal transition rates and ts/tv ratio = 2
    // Transversion rate = 1/ts_tv_ratio = 0.5 relative to transitions = 1
    let t = tn93(TN93Params {
      mu: 1.0,
      kappa1: 0.5, // transversions
      kappa2: 1.0, // C<->T same as A<->G
      pi: Some(pi),
      alphabet: AlphabetName::Nuc,
    })?;

    // Compare normalized Q matrices
    // Note: W is normalized by average rate, so relative rates should match
    let q_h = h.Q();
    let q_t = t.Q();
    assert_abs_diff_eq!(q_h, q_t, epsilon = 1e-14);

    // Also verify pi is the same
    assert_abs_diff_eq!(h.pi, t.pi, epsilon = 1e-14);

    Ok(())
  }

  /// TN93 is GTR with specific W matrix structure.
  /// GTR allows all 6 exchangeability parameters to vary independently.
  /// TN93 constrains: all transversions equal (kappa1), C<->T = kappa2, A<->G = 1.
  ///
  /// This test verifies the final link in the hierarchy: TN93 -> GTR.
  #[test]
  fn test_gtr_tn93_equals_gtr_with_structured_w() -> Result<(), Report> {
    let pi = array![0.1, 0.2, 0.3, 0.4];
    let kappa1 = 0.5; // transversion rate
    let kappa2 = 1.5; // C<->T rate (A<->G = 1 is reference)

    let t = tn93(TN93Params {
      mu: 1.0,
      kappa1,
      kappa2,
      pi: Some(pi.clone()),
      alphabet: AlphabetName::Nuc,
    })?;

    // GTR W matrix matching TN93 structure:
    // States: A=0, C=1, G=2, T=3
    // Transitions: A<->G (indices 0,2) = 1.0, C<->T (indices 1,3) = kappa2
    // Transversions: all others = kappa1
    #[rustfmt::skip]
    let w = array![
      [0.0,    kappa1, 1.0,    kappa1],  // A row: A-C=tv, A-G=ts(ref), A-T=tv
      [kappa1, 0.0,    kappa1, kappa2],  // C row: C-A=tv, C-G=tv, C-T=ts
      [1.0,    kappa1, 0.0,    kappa1],  // G row: G-A=ts(ref), G-C=tv, G-T=tv
      [kappa1, kappa2, kappa1, 0.0   ]   // T row: T-A=tv, T-C=ts, T-G=tv
    ];

    let alphabet = Alphabet::new(AlphabetName::Nuc).expect("Nuc alphabet should be valid");
    let n_states = alphabet.n_canonical();

    let g = GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: Some(w),
      pi,
    })?;

    // Compare normalized Q matrices
    let q_t = t.Q();
    let q_g = g.Q();
    assert_abs_diff_eq!(q_t, q_g, epsilon = 1e-14);

    // Also verify pi is the same
    assert_abs_diff_eq!(t.pi, g.pi, epsilon = 1e-14);

    Ok(())
  }
}
