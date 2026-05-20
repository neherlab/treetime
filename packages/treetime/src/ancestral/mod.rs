pub mod fitch;
pub mod fitch_indel;
pub mod fitch_sub;
pub(crate) mod gtr_inference;
#[allow(dead_code)]
pub(crate) mod gtr_inference_dense;
pub mod marginal;
pub mod params;

#[cfg(test)]
mod __tests__;
