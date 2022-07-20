//! # Least Squares
//!
//! Compute a least-squares solution to the equation Ax = b.
//! Compute a vector x such that the 2-norm ``|b - A x|`` is minimized.
//!
//! Finding the least squares solutions is implemented as traits, meaning
//! that to solve `A x = b` for a matrix `A` and a RHS `b`, we call
//! `let result = A.least_squares(&b);`. This returns a `result` of
//! type `LeastSquaresResult`, the solution for the least square problem
//! is in `result.solution`.
//!
//! There are three traits, `LeastSquaresSvd` with the method `least_squares`,
//! which operates on immutable references, `LeastSquaresInto` with the method
//! `least_squares_into`, which takes ownership over both the array `A` and the
//! RHS `b` and `LeastSquaresSvdInPlace` with the method `least_squares_in_place`,
//! which operates on mutable references for `A` and `b` and destroys these when
//! solving the least squares problem. `LeastSquaresSvdInto` and
//! `LeastSquaresSvdInPlace` avoid an extra allocation for `A` and `b` which
//! `LeastSquaresSvd` has do perform to preserve the values in `A` and `b`.
//!
//! All methods use the Lapacke family of methods `*gelsd` which solves the least
//! squares problem using the SVD with a divide-and-conquer strategy.
//!
//! The traits are implemented for value types `f32`, `f64`, `c32` and `c64`
//! and vector or matrix right-hand-sides (`ArrayBase<S, Ix1>` or `ArrayBase<S, Ix2>`).
//!
//! ## Example
//! ```rust
//! use approx::AbsDiffEq; // for abs_diff_eq
//! use ndarray::{array, Array1, Array2};
//! use ndarray_linalg::{LeastSquaresSvd, LeastSquaresSvdInto, LeastSquaresSvdInPlace};
//!
//! let a: Array2<f64> = array![
//!     [1., 1., 1.],
//!     [2., 3., 4.],
//!     [3., 5., 2.],
//!     [4., 2., 5.],
//!     [5., 4., 3.]
//! ];
//! // solving for a single right-hand side
//! let b: Array1<f64> = array![-10., 12., 14., 16., 18.];
//! let expected: Array1<f64> = array![2., 1., 1.];
//! let result = a.least_squares(&b).unwrap();
//! assert!(result.solution.abs_diff_eq(&expected, 1e-12));
//!
//! // solving for two right-hand sides at once
//! let b_2: Array2<f64> =
//!     array![[-10., -3.], [12., 14.], [14., 12.], [16., 16.], [18., 16.]];
//! let expected_2: Array2<f64> = array![[2., 1.], [1., 1.], [1., 2.]];
//! let result_2 = a.least_squares(&b_2).unwrap();
//! assert!(result_2.solution.abs_diff_eq(&expected_2, 1e-12));
//!
//! // using `least_squares_in_place` which overwrites its arguments
//! let mut a_3 = a.clone();
//! let mut b_3 = b.clone();
//! let result_3 = a_3.least_squares_in_place(&mut b_3).unwrap();
//!
//! // using `least_squares_into` which consumes its arguments
//! let result_4 = a.least_squares_into(b).unwrap();
//! // `a` and `b` have been moved, no longer valid
//! ```

use lax::*;
use ndarray::*;

use crate::error::*;
use crate::layout::*;
use crate::types::*;

/// Result of a LeastSquares computation
///
/// Takes two type parameters, `E`, the element type of the matrix
/// (one of `f32`, `f64`, `c32` or `c64`) and `I`, the dimension of
/// b in the equation `Ax = b` (one of `Ix1` or `Ix2`). If `I` is `Ix1`,
/// the  right-hand-side (RHS) is a `n x 1` column vector and the solution
/// is a `m x 1` column vector. If `I` is `Ix2`, the RHS is a `n x k` matrix
/// (which can be seen as solving `Ax = b` k times for different b) and
/// the solution is a `m x k` matrix.
#[derive(Debug, Clone)]
pub struct LeastSquaresResult<E: Scalar, I: Dimension> {
    /// The singular values of the matrix A in `Ax = b`
    pub singular_values: Array1<E::Real>,
    /// The solution vector or matrix `x` which is the best
    /// solution to `Ax = b`, i.e. minimizing the 2-norm `||b - Ax||`
    pub solution: Array<E, I>,
    /// The rank of the matrix A in `Ax = b`
    pub rank: i32,
    /// If n < m and rank(A) == n, the sum of squares
    /// If b is a (m x 1) vector, this is a 0-dimensional array (single value)
    /// If b is a (m x k) matrix, this is a (k x 1) column vector
    pub residual_sum_of_squares: Option<Array<E::Real, I::Smaller>>,
}
/// Solve least squares for immutable references
pub trait LeastSquaresSvd<D, E, I>
where
    D: Data<Elem = E>,
    E: Scalar + Lapack,
    I: Dimension,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(&rhs)`. `A` and `rhs`
    /// are unchanged.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares(&self, rhs: &ArrayBase<D, I>) -> Result<LeastSquaresResult<E, I>>;
}

/// Solve least squares for owned matrices
pub trait LeastSquaresSvdInto<D, E, I>
where
    D: Data<Elem = E>,
    E: Scalar + Lapack,
    I: Dimension,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(rhs)`, consuming both `A`
    /// and `rhs`. This uses the memory location of `A` and
    /// `rhs`, which avoids some extra memory allocations.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares_into(self, rhs: ArrayBase<D, I>) -> Result<LeastSquaresResult<E, I>>;
}

/// Solve least squares for mutable references, overwriting
/// the input fields in the process
pub trait LeastSquaresSvdInPlace<D, E, I>
where
    D: Data<Elem = E>,
    E: Scalar + Lapack,
    I: Dimension,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(&mut rhs)`, overwriting both `A`
    /// and `rhs`. This uses the memory location of `A` and
    /// `rhs`, which avoids some extra memory allocations.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares_in_place(
        &mut self,
        rhs: &mut ArrayBase<D, I>,
    ) -> Result<LeastSquaresResult<E, I>>;
}

/// Solve least squares for immutable references and a single
/// column vector as a right-hand side.
/// `E` is one of `f32`, `f64`, `c32`, `c64`. `D1`, `D2` can be any
/// valid representation for `ArrayBase` (over `E`).
impl<E, D1, D2> LeastSquaresSvd<D2, E, Ix1> for ArrayBase<D1, Ix2>
where
    E: Scalar + Lapack,
    D1: Data<Elem = E>,
    D2: Data<Elem = E>,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(&rhs)`, where `rhs` is a
    /// single column vector. `A` and `rhs` are unchanged.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares(&self, rhs: &ArrayBase<D2, Ix1>) -> Result<LeastSquaresResult<E, Ix1>> {
        let a = self.to_owned();
        let b = rhs.to_owned();
        a.least_squares_into(b)
    }
}

/// Solve least squares for immutable references and matrix
/// (=mulitipe vectors) as a right-hand side.
/// `E` is one of `f32`, `f64`, `c32`, `c64`. `D1`, `D2` can be any
/// valid representation for `ArrayBase` (over `E`).
impl<E, D1, D2> LeastSquaresSvd<D2, E, Ix2> for ArrayBase<D1, Ix2>
where
    E: Scalar + Lapack,
    D1: Data<Elem = E>,
    D2: Data<Elem = E>,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(&rhs)`, where `rhs` is
    /// matrix. `A` and `rhs` are unchanged.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares(&self, rhs: &ArrayBase<D2, Ix2>) -> Result<LeastSquaresResult<E, Ix2>> {
        let a = self.to_owned();
        let b = rhs.to_owned();
        a.least_squares_into(b)
    }
}

/// Solve least squares for owned values and a single
/// column vector as a right-hand side. The matrix and the RHS
/// vector are consumed.
///
/// `E` is one of `f32`, `f64`, `c32`, `c64`. `D` can be any
/// valid representation for `ArrayBase`.
impl<E, D1, D2> LeastSquaresSvdInto<D2, E, Ix1> for ArrayBase<D1, Ix2>
where
    E: Scalar + Lapack,
    D1: DataMut<Elem = E>,
    D2: DataMut<Elem = E>,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(rhs)`, where `rhs` is a
    /// single column vector. `A` and `rhs` are consumed.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares_into(
        mut self,
        mut rhs: ArrayBase<D2, Ix1>,
    ) -> Result<LeastSquaresResult<E, Ix1>> {
        self.least_squares_in_place(&mut rhs)
    }
}

/// Solve least squares for owned values and a matrix
/// as a right-hand side. The matrix and the RHS matrix
/// are consumed.
///
/// `E` is one of `f32`, `f64`, `c32`, `c64`. `D1`, `D2` can be any
/// valid representation for `ArrayBase` (over `E`).
impl<E, D1, D2> LeastSquaresSvdInto<D2, E, Ix2> for ArrayBase<D1, Ix2>
where
    E: Scalar + Lapack,
    D1: DataMut<Elem = E>,
    D2: DataMut<Elem = E>,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(rhs)`, where `rhs` is a
    /// matrix. `A` and `rhs` are consumed.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares_into(
        mut self,
        mut rhs: ArrayBase<D2, Ix2>,
    ) -> Result<LeastSquaresResult<E, Ix2>> {
        self.least_squares_in_place(&mut rhs)
    }
}

/// Solve least squares for mutable references and a vector
/// as a right-hand side. Both values are overwritten in the
/// call.
///
/// `E` is one of `f32`, `f64`, `c32`, `c64`. `D1`, `D2` can be any
/// valid representation for `ArrayBase` (over `E`).
impl<E, D1, D2> LeastSquaresSvdInPlace<D2, E, Ix1> for ArrayBase<D1, Ix2>
where
    E: Scalar + Lapack,
    D1: DataMut<Elem = E>,
    D2: DataMut<Elem = E>,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(rhs)`, where `rhs` is a
    /// vector. `A` and `rhs` are overwritten in the call.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares_in_place(
        &mut self,
        rhs: &mut ArrayBase<D2, Ix1>,
    ) -> Result<LeastSquaresResult<E, Ix1>> {
        if self.shape()[0] != rhs.shape()[0] {
            return Err(ShapeError::from_kind(ErrorKind::IncompatibleShape).into());
        }
        let (m, n) = (self.shape()[0], self.shape()[1]);
        if n > m {
            // we need a new rhs b/c it will be overwritten with the solution
            // for which we need `n` entries
            let mut new_rhs = Array1::<E>::zeros((n,));
            new_rhs.slice_mut(s![0..m]).assign(rhs);
            compute_least_squares_srhs(self, &mut new_rhs)
        } else {
            compute_least_squares_srhs(self, rhs)
        }
    }
}

fn compute_least_squares_srhs<E, D1, D2>(
    a: &mut ArrayBase<D1, Ix2>,
    rhs: &mut ArrayBase<D2, Ix1>,
) -> Result<LeastSquaresResult<E, Ix1>>
where
    E: Scalar + Lapack,
    D1: DataMut<Elem = E>,
    D2: DataMut<Elem = E>,
{
    let LeastSquaresOutput::<E> {
        singular_values,
        rank,
    } = E::least_squares(
        a.layout()?,
        a.as_allocated_mut()?,
        rhs.as_slice_memory_order_mut()
            .ok_or_else(|| LinalgError::MemoryNotCont)?,
    )?;

    let (m, n) = (a.shape()[0], a.shape()[1]);
    let solution = rhs.slice(s![0..n]).to_owned();
    let residual_sum_of_squares = compute_residual_scalar(m, n, rank, &rhs);
    Ok(LeastSquaresResult {
        solution,
        singular_values: Array::from_shape_vec((singular_values.len(),), singular_values)?,
        rank,
        residual_sum_of_squares,
    })
}

fn compute_residual_scalar<E: Scalar, D: Data<Elem = E>>(
    m: usize,
    n: usize,
    rank: i32,
    b: &ArrayBase<D, Ix1>,
) -> Option<Array<E::Real, Ix0>> {
    if m < n || n != rank as usize {
        return None;
    }
    let mut arr: Array<E::Real, Ix0> = Array::zeros(());
    arr[()] = b.slice(s![n..]).mapv(|x| x.powi(2).abs()).sum();
    Some(arr)
}

/// Solve least squares for mutable references and a matrix
/// as a right-hand side. Both values are overwritten in the
/// call.
///
/// `E` is one of `f32`, `f64`, `c32`, `c64`. `D1`, `D2` can be any
/// valid representation for `ArrayBase` (over `E`).
impl<E, D1, D2> LeastSquaresSvdInPlace<D2, E, Ix2> for ArrayBase<D1, Ix2>
where
    E: Scalar + Lapack + LeastSquaresSvdDivideConquer_,
    D1: DataMut<Elem = E>,
    D2: DataMut<Elem = E>,
{
    /// Solve a least squares problem of the form `Ax = rhs`
    /// by calling `A.least_squares(rhs)`, where `rhs` is a
    /// matrix. `A` and `rhs` are overwritten in the call.
    ///
    /// `A` and `rhs` must have the same layout, i.e. they must
    /// be both either row- or column-major format, otherwise a
    /// `IncompatibleShape` error is raised.
    fn least_squares_in_place(
        &mut self,
        rhs: &mut ArrayBase<D2, Ix2>,
    ) -> Result<LeastSquaresResult<E, Ix2>> {
        if self.shape()[0] != rhs.shape()[0] {
            return Err(ShapeError::from_kind(ErrorKind::IncompatibleShape).into());
        }
        let (m, n) = (self.shape()[0], self.shape()[1]);
        if n > m {
            // we need a new rhs b/c it will be overwritten with the solution
            // for which we need `n` entries
            let k = rhs.shape()[1];
            let mut new_rhs = match self.layout()? {
                MatrixLayout::C { .. } => Array2::<E>::zeros((n, k)),
                MatrixLayout::F { .. } => Array2::<E>::zeros((n, k).f()),
            };
            new_rhs.slice_mut(s![0..m, ..]).assign(rhs);
            compute_least_squares_nrhs(self, &mut new_rhs)
        } else {
            compute_least_squares_nrhs(self, rhs)
        }
    }
}

fn compute_least_squares_nrhs<E, D1, D2>(
    a: &mut ArrayBase<D1, Ix2>,
    rhs: &mut ArrayBase<D2, Ix2>,
) -> Result<LeastSquaresResult<E, Ix2>>
where
    E: Scalar + Lapack,
    D1: DataMut<Elem = E>,
    D2: DataMut<Elem = E>,
{
    let a_layout = a.layout()?;
    let rhs_layout = rhs.layout()?;
    let LeastSquaresOutput::<E> {
        singular_values,
        rank,
    } = E::least_squares_nrhs(
        a_layout,
        a.as_allocated_mut()?,
        rhs_layout,
        rhs.as_allocated_mut()?,
    )?;

    let solution: Array2<E> = rhs.slice(s![..a.shape()[1], ..]).to_owned();
    let singular_values = Array::from_shape_vec((singular_values.len(),), singular_values)?;
    let (m, n) = (a.shape()[0], a.shape()[1]);
    let residual_sum_of_squares = compute_residual_array1(m, n, rank, &rhs);
    Ok(LeastSquaresResult {
        solution,
        singular_values,
        rank,
        residual_sum_of_squares,
    })
}

fn compute_residual_array1<E: Scalar, D: Data<Elem = E>>(
    m: usize,
    n: usize,
    rank: i32,
    b: &ArrayBase<D, Ix2>,
) -> Option<Array1<E::Real>> {
    if m < n || n != rank as usize {
        return None;
    }
    Some(
        b.slice(s![n.., ..])
            .mapv(|x| x.powi(2).abs())
            .sum_axis(Axis(0)),
    )
}

#[cfg(test)]
mod tests {
    use crate::{error::LinalgError, *};
    use approx::AbsDiffEq;
    use ndarray::*;

    //
    // Test that the different least squares traits work as intended on the
    // different array types.
    //
    //               | least_squares | ls_into | ls_in_place |
    // --------------+---------------+---------+-------------+
    // Array         | yes           | yes     | yes         |
    // ArcArray      | yes           | no      | no          |
    // CowArray      | yes           | yes     | yes         |
    // ArrayView     | yes           | no      | no          |
    // ArrayViewMut  | yes           | no      | yes         |
    //

    fn assert_result<D1: Data<Elem = f64>, D2: Data<Elem = f64>>(
        a: &ArrayBase<D1, Ix2>,
        b: &ArrayBase<D2, Ix1>,
        res: &LeastSquaresResult<f64, Ix1>,
    ) {
        assert_eq!(res.rank, 2);
        let b_hat = a.dot(&res.solution);
        let rssq = (b - &b_hat).mapv(|x| x.powi(2)).sum();
        assert!(res.residual_sum_of_squares.as_ref().unwrap()[()].abs_diff_eq(&rssq, 1e-12));
        assert!(res
            .solution
            .abs_diff_eq(&array![-0.428571428571429, 0.85714285714285], 1e-12));
    }

    #[test]
    fn on_arc() {
        let a: ArcArray2<f64> = array![[1., 2.], [4., 5.], [3., 4.]].into_shared();
        let b: ArcArray1<f64> = array![1., 2., 3.].into_shared();
        let res = a.least_squares(&b).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn on_cow() {
        let a = CowArray::from(array![[1., 2.], [4., 5.], [3., 4.]]);
        let b = CowArray::from(array![1., 2., 3.]);
        let res = a.least_squares(&b).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn on_view() {
        let a: Array2<f64> = array![[1., 2.], [4., 5.], [3., 4.]];
        let b: Array1<f64> = array![1., 2., 3.];
        let av = a.view();
        let bv = b.view();
        let res = av.least_squares(&bv).unwrap();
        assert_result(&av, &bv, &res);
    }

    #[test]
    fn on_view_mut() {
        let mut a: Array2<f64> = array![[1., 2.], [4., 5.], [3., 4.]];
        let mut b: Array1<f64> = array![1., 2., 3.];
        let av = a.view_mut();
        let bv = b.view_mut();
        let res = av.least_squares(&bv).unwrap();
        assert_result(&av, &bv, &res);
    }

    #[test]
    fn on_cow_view() {
        let a = CowArray::from(array![[1., 2.], [4., 5.], [3., 4.]]);
        let b: Array1<f64> = array![1., 2., 3.];
        let bv = b.view();
        let res = a.least_squares(&bv).unwrap();
        assert_result(&a, &bv, &res);
    }

    #[test]
    fn into_on_owned() {
        let a: Array2<f64> = array![[1., 2.], [4., 5.], [3., 4.]];
        let b: Array1<f64> = array![1., 2., 3.];
        let ac = a.clone();
        let bc = b.clone();
        let res = ac.least_squares_into(bc).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn into_on_arc() {
        let a: ArcArray2<f64> = array![[1., 2.], [4., 5.], [3., 4.]].into_shared();
        let b: ArcArray1<f64> = array![1., 2., 3.].into_shared();
        let a2 = a.clone();
        let b2 = b.clone();
        let res = a2.least_squares_into(b2).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn into_on_cow() {
        let a = CowArray::from(array![[1., 2.], [4., 5.], [3., 4.]]);
        let b = CowArray::from(array![1., 2., 3.]);
        let a2 = a.clone();
        let b2 = b.clone();
        let res = a2.least_squares_into(b2).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn into_on_owned_cow() {
        let a: Array2<f64> = array![[1., 2.], [4., 5.], [3., 4.]];
        let b = CowArray::from(array![1., 2., 3.]);
        let ac = a.clone();
        let b2 = b.clone();
        let res = ac.least_squares_into(b2).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn in_place_on_owned() {
        let a = array![[1., 2.], [4., 5.], [3., 4.]];
        let b = array![1., 2., 3.];
        let mut a2 = a.clone();
        let mut b2 = b.clone();
        let res = a2.least_squares_in_place(&mut b2).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn in_place_on_cow() {
        let a = CowArray::from(array![[1., 2.], [4., 5.], [3., 4.]]);
        let b = CowArray::from(array![1., 2., 3.]);
        let mut a2 = a.clone();
        let mut b2 = b.clone();
        let res = a2.least_squares_in_place(&mut b2).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn in_place_on_mut_view() {
        let a = array![[1., 2.], [4., 5.], [3., 4.]];
        let b = array![1., 2., 3.];
        let mut a2 = a.clone();
        let mut b2 = b.clone();
        let av = &mut a2.view_mut();
        let bv = &mut b2.view_mut();
        let res = av.least_squares_in_place(bv).unwrap();
        assert_result(&a, &b, &res);
    }

    #[test]
    fn in_place_on_owned_cow() {
        let a = array![[1., 2.], [4., 5.], [3., 4.]];
        let b = CowArray::from(array![1., 2., 3.]);
        let mut a2 = a.clone();
        let mut b2 = b.clone();
        let res = a2.least_squares_in_place(&mut b2).unwrap();
        assert_result(&a, &b, &res);
    }

    //
    // Testing error cases
    //
    #[test]
    fn incompatible_shape_error_on_mismatching_num_rows() {
        let a: Array2<f64> = array![[1., 2.], [4., 5.], [3., 4.]];
        let b: Array1<f64> = array![1., 2.];
        match a.least_squares(&b) {
            Err(LinalgError::Shape(e)) if e.kind() == ErrorKind::IncompatibleShape => {}
            _ => panic!("Should be raise IncompatibleShape"),
        }
    }
}
