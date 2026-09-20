#[cfg(test)]
mod __tests__;

use ndarray::{Array2, Array3, Axis};

pub fn matmul_3d(lhs: &Array3<f64>, rhs: &Array3<f64>) -> Array3<f64> {
  let lhs_i1ka = lhs.view().insert_axis(Axis(1));
  let rhs_jka = rhs.view().permuted_axes([1, 0, 2]);
  let rhs_1jka = rhs_jka.insert_axis(Axis(0));
  (&lhs_i1ka * &rhs_1jka).sum_axis(Axis(2))
}

pub fn matvec_3d(mat: &Array3<f64>, vec: &Array2<f64>) -> Array2<f64> {
  (mat * &vec.view().insert_axis(Axis(0))).sum_axis(Axis(1))
}
