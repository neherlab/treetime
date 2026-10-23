#[cfg(test)]
mod tests {
  use crate::text::{Weight, text_width};
  use helpers::rendered_right_edge;
  use rstest::rstest;
  use treetime_utils::cast::f64_to_i32_saturating;

  #[rustfmt::skip]
  #[rstest]
  #[case::regular(Weight::Regular, "normal")]
  #[case::bold(   Weight::Bold,    "bold")]
  #[trace]
  fn test_text_width_matches_the_advance_that_the_svg_renderer_draws(
    #[case] weight: Weight,
    #[case] svg_weight: &str,
  ) {
    let width = text_width("HHHHHHHHHH", 40.0, weight).unwrap();
    let rendered = rendered_right_edge("HHHHHHHHHHH", svg_weight) - rendered_right_edge("H", svg_weight);
    let difference = (f64_to_i32_saturating(width.round()) - rendered).abs();
    assert!(difference <= 1, "measured {width} px, rendered {rendered} px");
  }

  mod helpers {
    use crate::convert::svg_to_png;
    use tiny_skia::Pixmap;
    use treetime_utils::cast::u32_to_usize;

    pub(super) fn rendered_right_edge(text: &str, weight: &str) -> i32 {
      let svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="800" height="80"><rect width="800" height="80" fill="white"/><text x="10" y="60" font-family="Lato" font-size="40" font-weight="{weight}">{text}</text></svg>"#
      );
      let png = Pixmap::decode_png(&svg_to_png(&svg).unwrap()).unwrap();
      let width = u32_to_usize(png.width());
      let right_edge = png
        .pixels()
        .iter()
        .enumerate()
        .filter(|(_, pixel)| pixel.red() < 128)
        .map(|(index, _)| index % width)
        .max()
        .unwrap();
      i32::try_from(right_edge).unwrap()
    }
  }
}
