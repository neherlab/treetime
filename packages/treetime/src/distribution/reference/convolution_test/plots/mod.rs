mod error_plots;
mod functions_and_convolution;
mod pointwise_plots;
mod spatial_plots;
mod utils;

pub use crate::distribution::reference::convolution_test::plots::error_plots::{
  plot_absolute_error, plot_error_histogram, plot_tolerance_metrics,
};
pub use crate::distribution::reference::convolution_test::plots::functions_and_convolution::plot_functions_and_convolution;
pub use crate::distribution::reference::convolution_test::plots::pointwise_plots::{
  plot_derivative_errors, plot_pointwise_error_profiles,
};
pub use crate::distribution::reference::convolution_test::plots::spatial_plots::plot_spatial_profiles;
pub use crate::distribution::reference::convolution_test::plots::utils::{
  combined_range, expand_range, tolerance_label,
};
