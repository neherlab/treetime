#[cfg(test)]
mod tests {
  use crate::text::{Weight, text_width};
  use helpers::{drawn_right_edge, rendered_right_edge};
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

  #[rustfmt::skip]
  #[rstest]
  #[case::regular(Weight::Regular)]
  #[case::bold(   Weight::Bold)]
  #[trace]
  fn test_text_with_a_character_that_plex_sans_lacks_draws_as_wide_as_measured(#[case] weight: Weight) {
    let text = "H\u{207b}\u{207b}\u{207b}\u{207b}\u{207b}\u{207b}\u{207b}\u{207b}\u{207b}";
    let width = text_width(text, 40.0, weight).unwrap();
    let drawn = drawn_right_edge(&format!("{text}H"), weight) - drawn_right_edge("H", weight);
    let difference = (f64_to_i32_saturating(width.round()) - drawn).abs();
    assert!(difference <= 1, "measured {width} px, drawn {drawn} px");
  }

  mod helpers {
    use crate::convert::svg_to_png;
    use crate::font::register_chart_fonts;
    use crate::text::{Weight, draw_text};
    use plotters::prelude::*;
    use plotters::style::text_anchor::{HPos, Pos, VPos};
    use tiny_skia::Pixmap;
    use treetime_utils::cast::u32_to_usize;

    pub(super) fn drawn_right_edge(text: &str, weight: Weight) -> i32 {
      register_chart_fonts().unwrap();
      let mut svg = String::new();
      {
        let area = SVGBackend::with_string(&mut svg, (800, 80)).into_drawing_area();
        area.fill(&WHITE).unwrap();
        draw_text(
          &area,
          text,
          (10.0, 20.0),
          40.0,
          weight,
          &BLACK,
          Pos::new(HPos::Left, VPos::Top),
        )
        .unwrap();
        area.present().unwrap();
      }
      right_edge(&svg)
    }

    pub(super) fn rendered_right_edge(text: &str, weight: &str) -> i32 {
      right_edge(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="800" height="80"><rect width="800" height="80" fill="white"/><text x="10" y="60" font-family="IBM Plex Sans" font-size="40" font-weight="{weight}">{text}</text></svg>"#
      ))
    }

    fn right_edge(svg: &str) -> i32 {
      let png = Pixmap::decode_png(&svg_to_png(svg).unwrap()).unwrap();
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
