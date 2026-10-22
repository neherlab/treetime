#![expect(
  unsafe_code,
  reason = "declares the OpenBLAS introspection functions, which the C library defines"
)]

use std::ffi::{c_char, c_int};

unsafe extern "C" {
  pub fn openblas_get_config() -> *const c_char;
  pub fn openblas_get_corename() -> *const c_char;
  pub fn openblas_get_parallel() -> c_int;
  pub fn openblas_get_num_threads() -> c_int;
  pub fn openblas_get_num_procs() -> c_int;
}
