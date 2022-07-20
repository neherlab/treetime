use ndarray::*;
use ndarray_linalg::*;

#[test]
fn diag_1d() {
    let d = arr1(&[1.0, 2.0]);
    let v = arr1(&[1.0, 1.0]);
    let dv = d.into_diagonal().apply(&v);
    assert_close_l2!(&dv, &arr1(&[1.0, 2.0]), 1e-7);
}

#[test]
fn diag_2d() {
    let d = arr1(&[1.0, 2.0]);
    let m = arr2(&[[1.0, 1.0], [1.0, 1.0]]);
    let dm = d.into_diagonal().apply2(&m);
    println!("dm = {:?}", dm);
    assert_close_l2!(&dm, &arr2(&[[1.0, 1.0], [2.0, 2.0]]), 1e-7);
}

#[test]
fn diag_2d_multi() {
    let d = arr1(&[1.0, 2.0]);
    let m = arr2(&[[1.0, 1.0], [1.0, 1.0]]);
    let dm = d.into_diagonal().apply2_into(m);
    println!("dm = {:?}", dm);
    assert_close_l2!(&dm, &arr2(&[[1.0, 1.0], [2.0, 2.0]]), 1e-7);
}
