#[cfg(test)]
mod tests {
  use crate::convert::{svg_to_pdf, svg_to_png};
  use eyre::Report;
  use helpers::{SVG, non_white_pixel_count, pdf_text};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tiny_skia::Pixmap;
  use treetime_utils::assert_error;

  #[test]
  fn test_png_has_the_pixel_size_of_the_svg() {
    let png = Pixmap::decode_png(&svg_to_png(SVG).unwrap()).unwrap();
    assert_eq!((300, 200), (png.width(), png.height()));
  }

  #[test]
  fn test_png_draws_the_text() {
    let png = Pixmap::decode_png(&svg_to_png(SVG).unwrap()).unwrap();
    assert!(non_white_pixel_count(&png) > 100);
  }

  #[test]
  fn test_pdf_page_is_three_quarters_of_a_point_per_pixel() {
    let pdf = pdf_text(&svg_to_pdf(SVG).unwrap());
    assert!(pdf.starts_with("%PDF-"));
    assert!(pdf.contains("/MediaBox[0 0 225 150]"), "{pdf}");
  }

  #[test]
  fn test_pdf_embeds_plex_sans_subsets() {
    let pdf = pdf_text(&svg_to_pdf(SVG).unwrap());
    assert!(pdf.contains("+IBMPlexSans"), "{pdf}");
    assert!(pdf.contains("+IBMPlexSans-Bold"), "{pdf}");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::png(svg_to_png)]
  #[case::pdf(svg_to_pdf)]
  fn test_conversion_gives_the_same_bytes_every_time(#[case] convert: fn(&str) -> Result<Vec<u8>, Report>) {
    assert_eq!(convert(SVG).unwrap(), convert(SVG).unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::png(svg_to_png)]
  #[case::pdf(svg_to_pdf)]
  fn test_conversion_of_invalid_svg_is_an_error(#[case] convert: fn(&str) -> Result<Vec<u8>, Report>) {
    assert_error!(
      convert("<svg"),
      "When reading the SVG chart: SVG data parsing failed cause the document does not have a root node"
    );
  }

  mod helpers {
    use indoc::indoc;
    use tiny_skia::Pixmap;

    pub(super) const SVG: &str = indoc! {r##"
      <svg width="300" height="200" viewBox="0 0 300 200" xmlns="http://www.w3.org/2000/svg">
      <rect x="0" y="0" width="300" height="200" fill="#FFFFFF"/>
      <text x="20" y="60" font-family="IBM Plex Sans, sans-serif" font-size="24">Clock regression</text>
      <text x="20" y="120" font-family="IBM Plex Sans, sans-serif" font-size="24" font-weight="bold">R² 0.98 χ²</text>
      </svg>
    "##};

    pub(super) fn non_white_pixel_count(pixmap: &Pixmap) -> usize {
      pixmap
        .pixels()
        .iter()
        .filter(|pixel| (pixel.red(), pixel.green(), pixel.blue()) != (255, 255, 255))
        .count()
    }

    pub(super) fn pdf_text(pdf: &[u8]) -> String {
      String::from_utf8_lossy(pdf).into_owned()
    }
  }
}
