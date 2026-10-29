use crate::font::{CHART_FONT_FAMILY, PLEX_SANS_BOLD, PLEX_SANS_REGULAR};
use ab_glyph::{Font, FontRef, ScaleFont};
use eyre::Report;
use num_traits::ToPrimitive;
use plotters::style::text_anchor::Pos;
use plotters::style::{Color, FontDesc, FontFamily, FontStyle, TextStyle};
use treetime_utils::make_internal_report;

const PLOTTERS_SIZE_PER_PIXEL: f64 = 1.24;

pub fn text_width(text: &str, size: f64, weight: Weight) -> Result<f64, Report> {
  let font = FontRef::try_from_slice(weight.font_bytes())
    .map_err(|_invalid| make_internal_report!("The embedded chart font {weight:?} is not a valid font"))?;
  let units_per_em = font
    .units_per_em()
    .ok_or_else(|| make_internal_report!("The embedded chart font {weight:?} has no units per em"))?;
  let size = size
    .to_f32()
    .ok_or_else(|| make_internal_report!("Text size {size} is outside the range of the font scale"))?;
  let scaled = font.as_scaled(size * font.height_unscaled() / units_per_em);
  let glyphs = text.chars().map(|c| scaled.glyph_id(c)).collect::<Vec<_>>();
  let advances: f32 = glyphs.iter().map(|&glyph| scaled.h_advance(glyph)).sum();
  let kerning: f32 = glyphs
    .iter()
    .zip(glyphs.iter().skip(1))
    .map(|(&left, &right)| scaled.kern(left, right))
    .sum();
  Ok(f64::from(advances + kerning))
}

pub fn text_style<C: Color>(size: f64, weight: Weight, color: &C, anchor: Pos) -> TextStyle<'static> {
  FontDesc::new(
    FontFamily::Name(CHART_FONT_FAMILY),
    size * PLOTTERS_SIZE_PER_PIXEL,
    weight.font_style(),
  )
  .color(color)
  .pos(anchor)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Weight {
  Regular,
  Bold,
}

impl Weight {
  fn font_bytes(self) -> &'static [u8] {
    match self {
      Self::Regular => PLEX_SANS_REGULAR,
      Self::Bold => PLEX_SANS_BOLD,
    }
  }

  fn font_style(self) -> FontStyle {
    match self {
      Self::Regular => FontStyle::Normal,
      Self::Bold => FontStyle::Bold,
    }
  }
}
