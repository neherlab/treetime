//! Rust binding to Intel-MKL including
//!
//! - [Vector Mathematical Functions (mkl_vml.h)](https://software.intel.com/en-us/onemkl-developer-reference-c-vector-mathematical-functions)
//! - [Statistical Functions (mkl_vsl.h)](https://software.intel.com/en-us/onemkl-developer-reference-c-statistical-functions)
//!
//! Other parts of Intel-MKL is served via
//!
//! - [blas-sys](https://crates.io/crates/blas-sys)
//! - [lapack-sys](https://crates.io/crates/lapack-sys)
//! - [lapacke-sys](https://crates.io/crates/lapacke-sys)
//! - [fftw-sys](https://crates.io/crates/fftw-sys)
//!
#![allow(
    improper_ctypes,
    non_upper_case_globals,
    non_camel_case_types,
    non_snake_case
)]

extern crate intel_mkl_src;

include!("mkl.rs");

// Test linking
#[cfg(test)]
mod tests {
    use super::*;
    use approx::ulps_eq;
    use rand::distributions::{Distribution, Uniform};
    use std::ffi::c_void;

    fn gen_rand_array(n: usize) -> Vec<f64> {
        let mut rng = rand::thread_rng();
        let between = Uniform::from(0.0..2.0 * std::f64::consts::PI);
        let mut buf = vec![0.0; n];
        for val in buf.iter_mut() {
            *val = between.sample(&mut rng);
        }
        buf
    }

    #[test]
    fn cos() {
        let n = 1024;
        let a = gen_rand_array(n);
        let mut b = vec![0.0_f64; n];
        unsafe {
            vdCos(n as i32, a.as_ptr(), b.as_mut_ptr());
        }
        for i in 0..n {
            ulps_eq!(b[i], a[i].cos(), max_ulps = 4, epsilon = std::f64::EPSILON);
        }
    }

    #[test]
    fn new_stream() {
        let mut stream: *mut c_void = std::ptr::null_mut();
        unsafe {
            vslNewStream(
                &mut stream as *mut *mut c_void,
                VSL_BRNG_MT19937 as i32,
                777,
            );
        }
    }
}
