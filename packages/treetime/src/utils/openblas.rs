#![allow(
  trivial_numeric_casts,
  unsafe_code,
  clippy::multiple_unsafe_ops_per_block,
  clippy::undocumented_unsafe_blocks
)]

use crate::io::json::{JsonPretty, json_write_str};
use serde::{Deserialize, Serialize};
use std::ffi::{CStr, c_char, c_int};

unsafe extern "C" {
  fn openblas_get_config() -> *const c_char;
  fn openblas_get_corename() -> *const c_char;
  fn openblas_get_parallel() -> c_int;
  fn openblas_get_num_threads() -> c_int;
  fn openblas_get_num_procs() -> c_int;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ParallelMode {
  Sequential,       // 0
  ParallelPlatform, // 1
  ParallelOpenMP,   // 2
  Unknown(i32),     // For unsupported or undefined values
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OpenBlasInfo {
  pub config: Option<String>,
  pub core_name: Option<String>,
  pub parallel_mode: ParallelMode,
  pub num_threads: i32,
  pub num_procs: i32,
}

pub fn print_openblas_info() {
  eprintln!("{}", get_openblas_info_str());
}

pub fn get_openblas_info_str() -> String {
  json_write_str(&get_openblas_info(), JsonPretty(true)).unwrap()
}

pub fn get_openblas_info() -> OpenBlasInfo {
  unsafe {
    let config = get_c_string(openblas_get_config());
    let corename = get_c_string(openblas_get_corename());
    let parallel_mode = match openblas_get_parallel() {
      0 => ParallelMode::Sequential,
      1 => ParallelMode::ParallelPlatform,
      2 => ParallelMode::ParallelOpenMP,
      other => ParallelMode::Unknown(other),
    };
    let num_threads = openblas_get_num_threads();
    let num_procs = openblas_get_num_procs();

    OpenBlasInfo {
      config,
      core_name: corename,
      parallel_mode,
      num_threads: num_threads as i32,
      num_procs: num_procs as i32,
    }
  }
}

fn get_c_string(ptr: *const c_char) -> Option<String> {
  (!ptr.is_null()).then(move || unsafe { CStr::from_ptr(ptr) }.to_string_lossy().into_owned())
}
