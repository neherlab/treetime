use ndarray::{Array2, Array3, Axis};

/// Batched matrix multiply over the third axis.
///
/// $C_{ija} = \sum_k A_{ika} \, B_{kja}$
///
/// Equivalent to `np.einsum('ika,kja->ija', A, B)`.
/// Each 2D slice `[:,:,a]` is an independent matrix multiply.
///
/// Implementation: expand `A` to `(i,1,k,a)`, permute `B` to `(j,k,a)` then
/// expand to `(1,j,k,a)`, broadcast-multiply to `(i,j,k,a)`, sum over `k`.
pub fn matmul_3d(lhs: &Array3<f64>, rhs: &Array3<f64>) -> Array3<f64> {
  let lhs_i1ka = lhs.view().insert_axis(Axis(1));
  let rhs_jka = rhs.view().permuted_axes([1, 0, 2]);
  let rhs_1jka = rhs_jka.insert_axis(Axis(0));
  (&lhs_i1ka * &rhs_1jka).sum_axis(Axis(2))
}

/// Batched matrix-vector multiply over the second axis.
///
/// $r_{ia} = \sum_j M_{ija} \, v_{ja}$
///
/// Equivalent to `np.einsum('ija,ja->ia', M, v)`.
/// Each 2D slice `M[:,:,a]` multiplies the corresponding column `v[:,a]`.
///
/// Implementation: expand `v` from `(j,a)` to `(1,j,a)`, broadcast-multiply
/// with `M` to `(i,j,a)`, sum over `j`.
pub fn matvec_3d(mat: &Array3<f64>, vec: &Array2<f64>) -> Array2<f64> {
  (mat * &vec.view().insert_axis(Axis(0))).sum_axis(Axis(1))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use ndarray::{Array1, Array2, Array3, s};
  use proptest::prelude::*;

  #[test]
  fn test_matmul_3d_matches_per_site_dot() {
    let a = Array3::from_shape_fn((2, 3, 2), |(i, j, a)| (i * 6 + j * 2 + a + 1) as f64);
    let b = Array3::from_shape_fn((3, 2, 2), |(k, j, a)| (k * 4 + j * 2 + a + 10) as f64);

    let result = matmul_3d(&a, &b);

    for site in 0..2 {
      let a_slice: Array2<f64> = a.slice(s![.., .., site]).to_owned();
      let b_slice: Array2<f64> = b.slice(s![.., .., site]).to_owned();
      let expected_slice = a_slice.dot(&b_slice);
      let result_slice: Array2<f64> = result.slice(s![.., .., site]).to_owned();
      pretty_assert_ulps_eq!(result_slice, expected_slice, max_ulps = 4);
    }
  }

  #[test]
  fn test_matmul_3d_nonsymmetric() {
    let a = Array3::from_shape_fn((2, 2, 1), |(i, j, _)| if i == 0 && j == 1 { 1.0 } else { 0.0 });
    let b = Array3::from_shape_fn((2, 2, 1), |(i, j, _)| if i == 0 && j == 0 { 1.0 } else { 0.0 });

    let ab = matmul_3d(&a, &b);
    let ba = matmul_3d(&b, &a);

    assert_ne!(ab, ba, "matmul_3d must not be commutative for non-symmetric matrices");
  }

  #[test]
  fn test_matvec_3d_matches_per_site_dot() {
    let mat = Array3::from_shape_fn((3, 3, 2), |(i, j, a)| (i * 6 + j * 2 + a + 1) as f64);
    let vec = Array2::from_shape_fn((3, 2), |(j, a)| (j * 2 + a + 1) as f64);

    let result = matvec_3d(&mat, &vec);

    for site in 0..2 {
      let m_slice: Array2<f64> = mat.slice(s![.., .., site]).to_owned();
      let v_slice: Array1<f64> = vec.slice(s![.., site]).to_owned();
      let expected_slice = m_slice.dot(&v_slice);
      let result_slice: Array1<f64> = result.slice(s![.., site]).to_owned();
      pretty_assert_ulps_eq!(result_slice, expected_slice, max_ulps = 4);
    }
  }

  #[test]
  fn test_matmul_3d_identity() {
    let a = Array3::from_shape_fn((3, 3, 2), |(i, j, a)| (i * 6 + j * 2 + a + 1) as f64);
    let identity = Array3::from_shape_fn((3, 3, 2), |(i, j, _)| if i == j { 1.0 } else { 0.0 });

    let ai = matmul_3d(&a, &identity);
    let ia = matmul_3d(&identity, &a);

    pretty_assert_abs_diff_eq!(ai, a, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(ia, a, epsilon = 1e-10);
  }

  fn arb_array3(rows: usize, cols: usize, sites: usize) -> impl Strategy<Value = Array3<f64>> {
    prop::collection::vec(-10.0_f64..10.0, rows * cols * sites)
      .prop_map(move |v| Array3::from_shape_vec((rows, cols, sites), v).unwrap())
  }

  fn arb_array2(rows: usize, sites: usize) -> impl Strategy<Value = Array2<f64>> {
    prop::collection::vec(-10.0_f64..10.0, rows * sites)
      .prop_map(move |v| Array2::from_shape_vec((rows, sites), v).unwrap())
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// $(AB)C = A(BC)$
    #[test]
    fn test_prop_matmul_3d_associative(
      a in arb_array3(3, 4, 2),
      b in arb_array3(4, 3, 2),
      c in arb_array3(3, 2, 2),
    ) {
      let ab_c = matmul_3d(&matmul_3d(&a, &b), &c);
      let a_bc = matmul_3d(&a, &matmul_3d(&b, &c));
      prop_assert!(
        approx::abs_diff_eq!(&ab_c, &a_bc, epsilon = 1e-8),
        "(AB)C != A(BC)"
      );
    }

    /// $AI = IA = A$
    #[test]
    fn test_prop_matmul_3d_identity(a in arb_array3(3, 3, 2)) {
      let identity = Array3::from_shape_fn((3, 3, 2), |(i, j, _)| if i == j { 1.0 } else { 0.0 });
      let ai = matmul_3d(&a, &identity);
      let ia = matmul_3d(&identity, &a);
      prop_assert!(approx::abs_diff_eq!(&ai, &a, epsilon = 1e-10), "AI != A");
      prop_assert!(approx::abs_diff_eq!(&ia, &a, epsilon = 1e-10), "IA != A");
    }

    /// $A(u + v) = Au + Av$
    #[test]
    fn test_prop_matmul_3d_distributive(
      a in arb_array3(3, 4, 2),
      b in arb_array3(4, 3, 2),
      c in arb_array3(4, 3, 2),
    ) {
      let a_bpc = matmul_3d(&a, &(&b + &c));
      let ab_plus_ac = &matmul_3d(&a, &b) + &matmul_3d(&a, &c);
      prop_assert!(
        approx::abs_diff_eq!(&a_bpc, &ab_plus_ac, epsilon = 1e-8),
        "A(B+C) != AB + AC"
      );
    }

    /// $M(av) = a(Mv)$
    #[test]
    fn test_prop_matvec_3d_scalar_linearity(
      mat in arb_array3(3, 3, 2),
      vec in arb_array2(3, 2),
      scalar in -10.0_f64..10.0,
    ) {
      let m_sv = matvec_3d(&mat, &(&vec * scalar));
      let s_mv = &matvec_3d(&mat, &vec) * scalar;
      prop_assert!(
        approx::abs_diff_eq!(&m_sv, &s_mv, epsilon = 1e-8),
        "M(av) != a(Mv)"
      );
    }

    /// $M(u + v) = Mu + Mv$
    #[test]
    fn test_prop_matvec_3d_additive(
      mat in arb_array3(3, 3, 2),
      u in arb_array2(3, 2),
      v in arb_array2(3, 2),
    ) {
      let m_upv = matvec_3d(&mat, &(&u + &v));
      let mu_plus_mv = &matvec_3d(&mat, &u) + &matvec_3d(&mat, &v);
      prop_assert!(
        approx::abs_diff_eq!(&m_upv, &mu_plus_mv, epsilon = 1e-8),
        "M(u+v) != Mu + Mv"
      );
    }

    /// $Iv = v$
    #[test]
    fn test_prop_matvec_3d_identity(vec in arb_array2(3, 2)) {
      let identity = Array3::from_shape_fn((3, 3, 2), |(i, j, _)| if i == j { 1.0 } else { 0.0 });
      let result = matvec_3d(&identity, &vec);
      prop_assert!(approx::abs_diff_eq!(&result, &vec, epsilon = 1e-10), "Iv != v");
    }

    /// $(AB)v = A(Bv)$
    #[test]
    fn test_prop_matvec_3d_matmul_consistency(
      a in arb_array3(3, 4, 2),
      b in arb_array3(4, 3, 2),
      v in arb_array2(3, 2),
    ) {
      let ab_v = matvec_3d(&matmul_3d(&a, &b), &v);
      let a_bv = matvec_3d(&a, &matvec_3d(&b, &v));
      prop_assert!(
        approx::abs_diff_eq!(&ab_v, &a_bv, epsilon = 1e-8),
        "(AB)v != A(Bv)"
      );
    }
  }
}
