use ndarray::*;
use ndarray_linalg::*;

#[test]
fn generalize() {
    let a: Array3<f64> = random((3, 2, 4).f());
    let ans = a.clone();
    let a: Array3<f64> = convert::generalize(a);
    assert_eq!(a, ans);
}
