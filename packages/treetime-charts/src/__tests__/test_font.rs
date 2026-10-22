#[cfg(test)]
mod tests {
  use crate::font::{CHART_FONT_FAMILY, register_chart_fonts};
  use helpers::text_svg;
  use indoc::indoc;
  use plotters::prelude::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_chart_text_names_lato_with_a_sans_serif_fallback() {
    let svg = text_svg((CHART_FONT_FAMILY, 20).into_font().style(FontStyle::Bold));
    assert_eq!(
      indoc! {r##"
        <svg width="200" height="50" viewBox="0 0 200 50" xmlns="http://www.w3.org/2000/svg">
        <text x="10" y="10" dy="0.76em" text-anchor="start" font-family="Lato, sans-serif" font-size="16.129032258064516" opacity="1" fill="#000000" font-weight="bold">
        Rate: 1.5e-3 χ²
        </text>
        </svg>
      "##},
      svg
    );
  }

  #[test]
  fn test_text_in_the_default_family_draws_with_lato() {
    let lato = text_svg((CHART_FONT_FAMILY, 20).into_font());
    let default = text_svg(("sans-serif", 20).into_font());
    assert_eq!(
      lato.replace(r#"font-family="Lato, sans-serif""#, r#"font-family="sans-serif""#),
      default
    );
  }

  mod helpers {
    use super::*;

    pub(super) fn text_svg(font: FontDesc<'_>) -> String {
      register_chart_fonts().unwrap();
      let mut text = String::new();
      {
        let area = SVGBackend::with_string(&mut text, (200, 50)).into_drawing_area();
        area.draw(&Text::new("Rate: 1.5e-3 χ²", (10, 10), font)).unwrap();
        area.present().unwrap();
      }
      text
    }
  }
}
