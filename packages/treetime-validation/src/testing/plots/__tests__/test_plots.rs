#[cfg(test)]
mod tests {
  use crate::algorithms::ConvolutionAlgorithm;
  use crate::testing::framework::results::TestResult;
  use crate::testing::plots::plots::validation_plots;
  use crate::testing::runners::convolution::ConvolutionRunner;
  use crate::testing::runners::runner::TestRunner;
  use crate::testing::test_suites::gaussian::{GaussianTestCase, GaussianTestSuite};
  use crate::testing::test_suites::test_suites::ConvolutionTestSuite;
  use eyre::Report;
  use helpers::{gaussian_result, render_plot};
  use itertools::Itertools;
  use plotters::prelude::{IntoDrawingArea, SVGBackend};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_charts::font::register_chart_fonts;

  #[rustfmt::skip]
  #[rstest]
  #[case::functions_convolution(   "functions_convolution.svg")]
  #[case::absolute_error(          "absolute_error.svg")]
  #[case::tolerance_metrics(       "tolerance_metrics.svg")]
  #[case::pointwise_error_profiles("pointwise_error_profiles.svg")]
  #[case::derivative_errors(       "derivative_errors.svg")]
  #[case::spatial_profiles(        "spatial_profiles.svg")]
  #[case::error_histogram(         "error_histogram.svg")]
  #[trace]
  fn test_plots_every_text_uses_the_chart_font(#[case] file_name: &str) -> Result<(), Report> {
    let svg = render_plot(file_name, &gaussian_result()?)?;

    let text_tags = svg
      .split("<text")
      .skip(1)
      .map(|rest| rest.split_once('>').map_or(rest, |(attributes, _)| attributes))
      .collect_vec();
    let tags_without_chart_font = text_tags
      .iter()
      .filter(|tag| !tag.contains(r#"font-family="Lato, sans-serif""#))
      .collect_vec();

    assert!(!text_tags.is_empty(), "the plot {file_name} has no text");
    assert_eq!(Vec::<&&str>::new(), tags_without_chart_font);
    Ok(())
  }

  mod helpers {
    use super::{
      ConvolutionAlgorithm, ConvolutionRunner, ConvolutionTestSuite, GaussianTestCase, GaussianTestSuite,
      IntoDrawingArea, Report, SVGBackend, TestResult, TestRunner, register_chart_fonts, validation_plots,
    };

    pub(super) fn gaussian_result() -> Result<TestResult<GaussianTestCase>, Report> {
      let suite = GaussianTestSuite;
      let test_case = suite
        .create_test_cases()
        .into_iter()
        .next()
        .expect("the Gaussian suite has test cases");
      ConvolutionRunner::<GaussianTestSuite>::run_single_test(&suite, &test_case, ConvolutionAlgorithm::Direct)
    }

    pub(super) fn render_plot(file_name: &str, result: &TestResult<GaussianTestCase>) -> Result<String, Report> {
      register_chart_fonts()?;
      let plot = validation_plots()
        .into_iter()
        .find(|plot| plot.file_name == file_name)
        .expect("the plot is listed");
      let mut svg = String::new();
      {
        let root = SVGBackend::with_string(&mut svg, plot.size).into_drawing_area();
        (plot.draw)(result, &root)?;
        root.present()?;
      }
      Ok(svg)
    }
  }
}
