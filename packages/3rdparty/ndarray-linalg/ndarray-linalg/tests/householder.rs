use ndarray::*;
use ndarray_linalg::{krylov::*, *};

fn over<A: Scalar + Lapack>(rtol: A::Real) {
    const N: usize = 4;
    let a: Array2<A> = random((N, N * 2));

    // Terminate
    let (q, r) = householder(a.axis_iter(Axis(1)), N, rtol, Strategy::Terminate);
    let a_sub = a.slice(s![.., 0..N]);
    let qc: Array2<A> = conjugate(&q);
    assert_close_l2!(&qc.dot(&q), &Array::eye(N), rtol; "Check Q^H Q = I");
    assert_close_l2!(&q.dot(&r), &a_sub, rtol; "Check A = QR");

    // Skip
    let (q, r) = householder(a.axis_iter(Axis(1)), N, rtol, Strategy::Skip);
    let a_sub = a.slice(s![.., 0..N]);
    let qc: Array2<A> = conjugate(&q);
    assert_close_l2!(&qc.dot(&q), &Array::eye(N), rtol);
    assert_close_l2!(&q.dot(&r), &a_sub, rtol);

    // Full
    let (q, r) = householder(a.axis_iter(Axis(1)), N, rtol, Strategy::Full);
    let qc: Array2<A> = conjugate(&q);
    assert_close_l2!(&qc.dot(&q), &Array::eye(N), rtol);
    assert_close_l2!(&q.dot(&r), &a, rtol);
}

#[test]
fn over_f32() {
    over::<f32>(1e-5);
}
#[test]
fn over_f64() {
    over::<f64>(1e-9);
}
#[test]
fn over_c32() {
    over::<c32>(1e-5);
}
#[test]
fn over_c64() {
    over::<c64>(1e-9);
}

fn full<A: Scalar + Lapack>(rtol: A::Real) {
    const N: usize = 5;
    let a: Array2<A> = random((N, N));
    let (q, r) = householder(a.axis_iter(Axis(1)), N, rtol, Strategy::Terminate);
    let qc: Array2<A> = conjugate(&q);
    assert_close_l2!(&qc.dot(&q), &Array::eye(N), rtol; "Check Q^H Q = I");
    assert_close_l2!(&q.dot(&r), &a, rtol; "Check A = QR");
}

#[test]
fn full_f32() {
    full::<f32>(1e-5);
}
#[test]
fn full_f64() {
    full::<f64>(1e-9);
}
#[test]
fn full_c32() {
    full::<c32>(1e-5);
}
#[test]
fn full_c64() {
    full::<c64>(1e-9);
}

fn half<A: Scalar + Lapack>(rtol: A::Real) {
    const N: usize = 4;
    let a: Array2<A> = random((N, N / 2));
    let (q, r) = householder(a.axis_iter(Axis(1)), N, rtol, Strategy::Terminate);
    let qc: Array2<A> = conjugate(&q);
    assert_close_l2!(&qc.dot(&q), &Array::eye(N / 2), rtol; "Check Q^H Q = I");
    assert_close_l2!(&q.dot(&r), &a, rtol; "Check A = QR");
}

#[test]
fn half_f32() {
    half::<f32>(1e-5);
}
#[test]
fn half_f64() {
    half::<f64>(1e-9);
}
#[test]
fn half_c32() {
    half::<c32>(1e-5);
}
#[test]
fn half_c64() {
    half::<c64>(1e-9);
}
