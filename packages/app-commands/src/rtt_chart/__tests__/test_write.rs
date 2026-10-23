#[cfg(test)]
mod tests {
  use crate::rtt_chart::__tests__::test_support::tests::{estimated_model, tip};
  use crate::rtt_chart::write::{encode_chart, rtt_chart_svg, write_rtt_chart};
  use ab_glyph::{Font, FontRef};
  use app_output::output_format::ChartFormat;
  use helpers::{chart_svg, non_white_pixel_count, read_decompressed, results, svg_attribute_values, svg_texts};
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::fs;
  use tempfile::TempDir;
  use tiny_skia::Pixmap;
  use treetime_charts::font::{LATO_BOLD, LATO_REGULAR};
  use treetime_utils::assert_error;
  use treetime_utils::fmt::float::float_to_power_of_ten;

  #[test]
  fn test_rtt_chart_without_dated_samples_is_an_error() {
    assert_error!(
      rtt_chart_svg(
        &[tip(None, 0.01, false), tip(None, 0.02, true)],
        &estimated_model(0.002, -4.0, 0.99)
      ),
      "the root-to-tip chart needs at least one sample with a date"
    );
  }

  #[test]
  fn test_rtt_chart_svg_writes_every_text_in_lato() {
    let families: BTreeSet<String> = svg_attribute_values(&chart_svg(&results()), "font-family")
      .into_iter()
      .collect();
    assert_eq!(btreeset! {"Lato, sans-serif".to_owned()}, families);
  }

  #[test]
  fn test_rtt_chart_svg_writes_whole_pixel_font_sizes() {
    let sizes: BTreeSet<String> = svg_attribute_values(&chart_svg(&results()), "font-size")
      .into_iter()
      .collect();
    assert_eq!(
      btreeset! {"12".to_owned(), "14".to_owned(), "24".to_owned(), "28".to_owned()},
      sizes
    );
  }

  #[test]
  fn test_rtt_chart_svg_names_itself_for_assistive_technology() {
    let svg = chart_svg(&results());
    let title = svg
      .split_once("<title>")
      .and_then(|(_, rest)| rest.split_once("</title>"))
      .unwrap()
      .0;
    let description = svg
      .split_once("<desc>")
      .and_then(|(_, rest)| rest.split_once("</desc>"))
      .unwrap()
      .0;
    assert_eq!(
      (
        "Root-to-tip regression",
        "Clock rate 2.00 × 10⁻³ (substitutions per site per year). Root date 2000.0. R² 0.980. 4 samples in the \
         fit, 1 excluded by the clock filter",
        true
      ),
      (title, description, svg.contains(r#"<svg role="img""#))
    );
  }

  #[test]
  fn test_rtt_chart_svg_shows_the_statistics_legend_and_axis_titles() {
    let texts = svg_texts(&chart_svg(&results()));
    let expected = [
      "Root-to-tip regression",
      "4 samples in the fit, 1 excluded by the clock filter",
      "Clock rate",
      "2.00 × 10⁻³",
      "substitutions per site per year",
      "Root date",
      "2000.0",
      "R²",
      "0.980",
      "Samples in the fit",
      "Excluded by the clock filter",
      "Least-squares fit",
      "Root-to-tip divergence (substitutions per site)",
      "Sampling date",
      "Root 2000.0",
    ];
    let missing = expected
      .iter()
      .filter(|text| !texts.iter().any(|t| t == *text))
      .collect_vec();
    assert_eq!(Vec::<&&str>::new(), missing, "{texts:?}");
  }

  #[test]
  fn test_rtt_chart_svg_without_outliers_lists_no_excluded_samples() {
    let fitted_only = results().into_iter().filter(|result| !result.is_outlier).collect_vec();
    let texts = svg_texts(&chart_svg(&fitted_only));
    let excluded = texts.iter().filter(|text| text.contains("xcluded")).collect_vec();
    assert_eq!(Vec::<&String>::new(), excluded);
  }

  #[test]
  fn test_rtt_chart_png_is_1200_by_800_with_drawn_pixels() {
    let png = Pixmap::decode_png(&encode_chart(ChartFormat::Png, chart_svg(&results())).unwrap()).unwrap();
    assert_eq!((1200, 800), (png.width(), png.height()));
    assert!(non_white_pixel_count(&png) > 1000);
  }

  #[test]
  fn test_rtt_chart_pdf_is_a_900_by_600_page_with_a_lato_subset() {
    let pdf = encode_chart(ChartFormat::Pdf, chart_svg(&results())).unwrap();
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
      encode_chart(format, chart_svg(&results())).unwrap(),
      encode_chart(format, chart_svg(&results())).unwrap()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::svg("chart.svg",    ChartFormat::Svg)]
  #[case::png("chart.PNG",    ChartFormat::Png)]
  #[case::pdf("chart.pdf",    ChartFormat::Pdf)]
  #[case::gz( "chart.svg.gz", ChartFormat::Svg)]
  #[trace]
  fn test_rtt_chart_writes_the_format_of_the_file_extension(#[case] name: &str, #[case] format: ChartFormat) {
    let dir = TempDir::new().unwrap();
    let path = dir.path().join(name);

    write_rtt_chart(&results(), &estimated_model(0.002, -4.0, 0.99), &path).unwrap();

    let files = fs::read_dir(dir.path())
      .unwrap()
      .map(|entry| entry.unwrap().file_name().into_string().unwrap())
      .collect_vec();
    assert_eq!(
      (vec![name.to_owned()], encode_chart(format, chart_svg(&results())).unwrap()),
      (files, read_decompressed(&path))
    );
  }

  #[test]
  fn test_rtt_chart_characters_have_glyphs_in_both_fonts() {
    let numbers = [
      0.0,
      -0.5,
      1.0,
      2005.25,
      -2005.25,
      0.00123,
      -1.5e-12,
      6.02e23,
      f64::NAN,
      f64::INFINITY,
      f64::NEG_INFINITY,
    ];
    let texts = numbers
      .iter()
      .map(|&x| float_to_power_of_ten(x, 3))
      .chain(["\u{2212}12", "± 0.31", "\u{2190} Root date", "Root date \u{2192}"].map(str::to_owned))
      .chain(svg_texts(&chart_svg(&results())))
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
    use std::io::Read;
    use std::path::Path;
    use treetime::clock::rtt::ClockRegressionResult;
    use treetime_utils::io::file::read_file_with;

    pub(super) fn results() -> Vec<ClockRegressionResult> {
      [
        (2000.0, 0.001, false),
        (2003.5, 0.008, false),
        (2007.0, 0.017, false),
        (2010.0, 0.019, false),
        (2008.0, 0.03, true),
      ]
      .into_iter()
      .map(|(date, div, is_outlier)| tip(Some(date), div, is_outlier))
      .collect()
    }

    pub(super) fn chart_svg(results: &[ClockRegressionResult]) -> String {
      rtt_chart_svg(results, &estimated_model(0.002, -4.0, 0.99)).unwrap()
    }

    pub(super) fn read_decompressed(path: &Path) -> Vec<u8> {
      read_file_with(path, |mut reader| {
        let mut bytes = vec![];
        reader.read_to_end(&mut bytes)?;
        Ok(bytes)
      })
      .unwrap()
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
