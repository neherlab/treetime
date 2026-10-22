use crate::testing::framework::results::{TestResult, TestRunOutcome};
use crate::testing::framework::test_case::TestCase;
use crate::testing::plots::error_plots::{plot_absolute_error, plot_error_histogram, plot_tolerance_metrics};
use crate::testing::plots::functions_and_convolution::plot_functions_and_convolution;
use crate::testing::plots::pointwise_plots::{plot_derivative_errors, plot_pointwise_error_profiles};
use crate::testing::plots::spatial_plots::plot_spatial_profiles;
use eyre::{Report, WrapErr};
use plotters::coord::Shift;
use plotters::prelude::{DrawingArea, IntoDrawingArea, SVGBackend};
use std::fs;
use treetime_charts::font::register_chart_fonts;

pub(crate) fn generate_plot_outputs<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>]) -> Result<(), Report>
where
  T: TestCase,
{
  register_chart_fonts()?;
  for outcome in outcomes {
    if let TestRunOutcome::Success(result) = outcome {
      let result = result.as_ref();
      let algorithm_dir = format!("{}/{}/{}", output_dir, result.test_case.name(), result.algorithm);
      fs::create_dir_all(&algorithm_dir).wrap_err_with(|| format!("When creating directory '{algorithm_dir}'"))?;

      for plot in validation_plots() {
        let output_path = format!("{algorithm_dir}/{}", plot.file_name);
        let root = SVGBackend::new(&output_path, plot.size).into_drawing_area();
        (plot.draw)(result, &root)?;
        root
          .present()
          .wrap_err_with(|| format!("When writing plot '{output_path}'"))?;
      }
    }
  }
  Ok(())
}

pub(crate) fn validation_plots<T>() -> [ValidationPlot<T>; 7]
where
  T: TestCase,
{
  [
    ValidationPlot {
      file_name: "functions_convolution.svg",
      size: (960, 640),
      draw: plot_functions_and_convolution,
    },
    ValidationPlot {
      file_name: "absolute_error.svg",
      size: (960, 480),
      draw: plot_absolute_error,
    },
    ValidationPlot {
      file_name: "tolerance_metrics.svg",
      size: (800, 600),
      draw: plot_tolerance_metrics,
    },
    ValidationPlot {
      file_name: "pointwise_error_profiles.svg",
      size: (1200, 900),
      draw: plot_pointwise_error_profiles,
    },
    ValidationPlot {
      file_name: "derivative_errors.svg",
      size: (1200, 450),
      draw: plot_derivative_errors,
    },
    ValidationPlot {
      file_name: "spatial_profiles.svg",
      size: (1200, 900),
      draw: plot_spatial_profiles,
    },
    ValidationPlot {
      file_name: "error_histogram.svg",
      size: (960, 640),
      draw: plot_error_histogram,
    },
  ]
}

pub(crate) struct ValidationPlot<T: TestCase> {
  pub(crate) file_name: &'static str,
  pub(crate) size: (u32, u32),
  pub(crate) draw: fn(&TestResult<T>, &DrawingArea<SVGBackend<'_>, Shift>) -> Result<(), Report>,
}
