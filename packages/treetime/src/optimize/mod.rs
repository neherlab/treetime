#[cfg(test)]
mod __tests__;

pub(super) mod dense_eval;
pub mod dispatch;
pub(super) mod eval;
pub mod indel;
pub mod iteration;
pub mod likelihood;
pub(super) mod method_brent;
pub(super) mod method_newton;
pub mod params;
pub mod pipeline;
pub(super) mod run_loop;
pub(super) mod sparse_eval;
pub mod topology;
pub(super) mod zero_boundary;
