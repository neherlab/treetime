#[cfg(test)]
mod __tests__;

pub mod algo;
pub mod discrete_states;
pub mod fitch;
pub mod fitch_config;
pub mod likelihood;
pub mod marginal_core;
pub mod marginal_dense;
pub mod marginal_discrete;
pub(crate) mod marginal_helpers;
mod marginal_passes;
pub mod marginal_sparse;
pub mod optimization_contribution;
pub mod optimize_dense;
pub mod optimize_sparse;
pub mod payload;
pub mod timetree;
pub mod traits;
