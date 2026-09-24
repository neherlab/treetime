pub(crate) mod distribution_core;
pub(crate) mod distribution_ops;
pub(crate) mod policy;

pub use distribution_core::distribution::Distribution;
pub use distribution_core::formula::DistributionFormula;
pub use distribution_core::function::DistributionFunction;
pub use distribution_core::point::DistributionPoint;
pub use distribution_core::range::DistributionRange;
pub use distribution_ops::divide::distribution_division;
pub use distribution_ops::edge_convolution::convolve_across_edge;
pub use distribution_ops::mass_domain::rewindow_to_mass;
pub use distribution_ops::multiply::distribution_multiplication;
pub use distribution_ops::multiply_by_fn::distribution_multiply_by_fn;
pub use distribution_ops::product::distribution_product;
pub use policy::{NegLog, Plain, PolicyMarker, YAxisPolicy};
pub use treetime_grid::BoundaryBehavior;

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor(unsafe)]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
