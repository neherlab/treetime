#[cfg(test)]
mod tests {
  use crate::rtt_chart::{clock_regression_chart_svg, gather_points, write_clock_regression_charts};
  use ab_glyph::{Font, FontRef};
  use helpers::{chart_svg, clock_model, non_white_pixel_count, row, svg_attribute_values, svg_texts};
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::fs;
  use tempfile::TempDir;
  use tiny_skia::Pixmap;
  use treetime_charts::convert::ChartFormat;
  use treetime_charts::font::{LATO_BOLD, LATO_REGULAR};
  use treetime_utils::assert_error;
  use treetime_utils::fmt::float::float_to_significant_digits;

  #[test]
  fn test_rtt_chart_without_dated_samples_is_an_error() {
    assert_error!(
      gather_points(&[row(None, false), row(None, true)], &clock_model()),
      "the root-to-tip chart needs at least one sample with a date"
    );
  }

  #[test]
  fn test_rtt_chart_line_spans_the_dated_samples() {
    let points = gather_points(
      &[row(Some(2000.0), false), row(Some(2010.0), true), row(None, false)],
      &clock_model(),
    )
    .unwrap();
    assert_eq!(
      ((2000.0, 0.0), (2010.0, 0.02), 1, 1),
      (
        points.line[0],
        points.line[1],
        points.norm_points.len(),
        points.outlier_points.len()
      )
    );
  }

  #[test]
  fn test_rtt_chart_svg_writes_every_text_in_lato() {
    let families: BTreeSet<String> = svg_attribute_values(&chart_svg(), "font-family").into_iter().collect();
    assert_eq!(btreeset! {"Lato, sans-serif".to_owned()}, families);
  }

  #[test]
  fn test_rtt_chart_svg_has_the_title_axis_descriptions_and_legend() {
    let labels = svg_texts(&chart_svg())
      .into_iter()
      .filter(|text| text.parse::<f64>().is_err())
      .collect_vec();
    assert_eq!(
      vec![
        "Clock regression",
        "Div",
        "Date",
        "Samples (norm)",
        "Samples (outliers)",
        "Clock regression: div = 0.002t - 4",
        "tMRCA: 2000",
        "Rate: 0.002",
        "Intercept: -4",
        "R: 0.99",
        "R²: 0.9801",
        "χ²: 0.000e0",
      ],
      labels
    );
  }

  #[test]
  fn test_rtt_chart_png_is_1200_by_800_with_drawn_pixels() {
    let png = Pixmap::decode_png(&ChartFormat::Png.encode(&chart_svg()).unwrap()).unwrap();
    assert_eq!((1200, 800), (png.width(), png.height()));
    assert!(non_white_pixel_count(&png) > 1000);
  }

  #[test]
  fn test_rtt_chart_pdf_is_a_900_by_600_page_with_a_lato_subset() {
    let pdf = ChartFormat::Pdf.encode(&chart_svg()).unwrap();
    let pdf = String::from_utf8_lossy(&pdf);
    assert!(pdf.starts_with("%PDF-"));
    assert!(pdf.contains("/MediaBox[0 0 900 600]"), "{pdf}");
    assert!(pdf.contains("+Lato-Regular"), "{pdf}");
    assert!(pdf.contains("+Lato-Bold"), "{pdf}");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::svg(ChartFormat::Svg)]
  #[case::png(ChartFormat::Png)]
  #[case::pdf(ChartFormat::Pdf)]
  #[trace]
  fn test_rtt_chart_renders_the_same_bytes_twice(#[case] format: ChartFormat) {
    assert_eq!(
      format.encode(&chart_svg()).unwrap(),
      format.encode(&chart_svg()).unwrap()
    );
  }

  #[test]
  fn test_rtt_chart_writes_only_the_requested_formats() {
    let dir = TempDir::new().unwrap();
    let svg_path = dir.path().join("chart.svg");
    let pdf_path = dir.path().join("chart.pdf");
    let results = helpers::results();

    write_clock_regression_charts(
      &results,
      &clock_model(),
      &[
        (ChartFormat::Svg, svg_path.as_path()),
        (ChartFormat::Pdf, pdf_path.as_path()),
      ],
    )
    .unwrap();

    let svg = chart_svg();
    let files = fs::read_dir(dir.path())
      .unwrap()
      .map(|entry| entry.unwrap().file_name().into_string().unwrap())
      .sorted()
      .collect_vec();
    assert_eq!(
      (
        vec!["chart.pdf".to_owned(), "chart.svg".to_owned()],
        svg.clone().into_bytes(),
        ChartFormat::Pdf.encode(&svg).unwrap()
      ),
      (files, fs::read(&svg_path).unwrap(), fs::read(&pdf_path).unwrap())
    );
  }

  #[test]
  fn test_rtt_chart_characters_have_glyphs_in_both_fonts() {
    let numbers = [
      0.0,
      -0.0,
      0.5,
      -0.5,
      1.0,
      -1.0,
      2005.25,
      -2005.25,
      0.00123,
      -0.00123,
      1.5e-12,
      -1.5e-12,
      6.02e23,
      -6.02e23,
      f64::NAN,
      f64::INFINITY,
      f64::NEG_INFINITY,
    ];
    let texts = numbers
      .iter()
      .flat_map(|&x| (1..=6).map(move |digits| float_to_significant_digits(x, digits)))
      .chain(numbers.iter().map(|x| format!("{x:.3e}")))
      .chain(svg_texts(&chart_svg()))
      .collect_vec();
    let characters: BTreeSet<char> = texts.iter().flat_map(|text| text.chars()).collect();

    let missing = [("Lato Regular", LATO_REGULAR), ("Lato Bold", LATO_BOLD)]
      .into_iter()
      .flat_map(|(name, bytes)| {
        let font = FontRef::try_from_slice(bytes).unwrap();
        characters
          .iter()
          .filter(move |&&c| font.glyph_id(c).0 == 0)
          .map(move |&c| (name, c))
          .collect_vec()
      })
      .collect_vec();
    assert_eq!(Vec::<(&str, char)>::new(), missing, "characters: {characters:?}");
  }

  mod helpers {
    use super::*;
    use indoc::indoc;
    use treetime::clock::clock_model::ClockModel;
    use treetime::clock::rtt::ClockRegressionResult;
    use treetime_utils::io::json::json_read_str;

    pub(super) fn row(date: Option<f64>, is_outlier: bool) -> ClockRegressionResult {
      ClockRegressionResult {
        name: String::new(),
        div: 0.01,
        date,
        predicted_date: 2005.0,
        clock_deviation: None,
        is_outlier,
        is_leaf: true,
        date_source: None,
      }
    }

    pub(super) fn clock_model() -> ClockModel {
      json_read_str(indoc! {r#"{
        "clock_rate": 0.002,
        "intercept": -4.0,
        "stats": {
          "estimated": {
            "chisq": 0.0,
            "r_val": 0.99,
            "hessian": [[0.0, 0.0], [0.0, 0.0]],
            "cov": [[1e-8, 0.0], [0.0, 0.5]]
          }
        }
      }"#})
      .unwrap()
    }

    pub(super) fn results() -> Vec<ClockRegressionResult> {
      [
        (2000.0, 0.001, false),
        (2003.5, 0.008, false),
        (2007.0, 0.017, false),
        (2010.0, 0.019, false),
        (2008.0, 0.03, true),
      ]
      .into_iter()
      .map(|(date, div, is_outlier)| ClockRegressionResult {
        div,
        ..row(Some(date), is_outlier)
      })
      .collect()
    }

    pub(super) fn chart_svg() -> String {
      clock_regression_chart_svg(&results(), &clock_model()).unwrap()
    }

    pub(super) fn svg_texts(svg: &str) -> Vec<String> {
      svg
        .split("<text")
        .skip(1)
        .filter_map(|element| {
          let (_, rest) = element.split_once('>')?;
          let (content, _) = rest.split_once("</text>")?;
          Some(content.trim().to_owned())
        })
        .collect()
    }

    pub(super) fn svg_attribute_values(svg: &str, attribute: &str) -> Vec<String> {
      let prefix = format!(r#" {attribute}=""#);
      svg
        .split(prefix.as_str())
        .skip(1)
        .filter_map(|rest| rest.split_once('"').map(|(value, _)| value.to_owned()))
        .collect()
    }

    pub(super) fn non_white_pixel_count(pixmap: &Pixmap) -> usize {
      pixmap
        .pixels()
        .iter()
        .filter(|pixel| (pixel.red(), pixel.green(), pixel.blue()) != (255, 255, 255))
        .count()
    }
  }
}
