#[cfg(test)]
mod __tests__;

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
