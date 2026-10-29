use crate::font::{
  CHART_FONT_FAMILY, FALLBACK_BOLD, FALLBACK_FONT_FAMILY, FALLBACK_REGULAR, PLEX_SANS_BOLD, PLEX_SANS_REGULAR,
};
use ab_glyph::{Font, FontRef, GlyphId, PxScaleFont, ScaleFont};
use eyre::Report;
use itertools::Itertools;
use num_traits::ToPrimitive;
use plotters::coord::Shift;
use plotters::element::Text;
use plotters::prelude::{DrawingArea, DrawingBackend};
use plotters::style::text_anchor::{HPos, Pos};
use plotters::style::{Color, FontDesc, FontFamily, FontStyle, TextStyle};
use treetime_utils::cast::f64_to_i32_saturating;
use treetime_utils::make_internal_report;

const PLOTTERS_SIZE_PER_PIXEL: f64 = 1.24;

pub fn text_width(text: &str, size: f64, weight: Weight) -> Result<f64, Report> {
  Ok(text_runs(text, size, weight)?.iter().map(|run| run.width).sum())
}

pub fn draw_text<DB, C>(
  area: &DrawingArea<DB, Shift>,
  text: &str,
  (x, y): (f64, f64),
  size: f64,
  weight: Weight,
  color: &C,
  anchor: Pos,
) -> Result<(), Report>
where
  DB: DrawingBackend,
  DB::ErrorType: 'static,
  C: Color,
{
  let runs = text_runs(text, size, weight)?;
  if let [run] = runs.as_slice()
    && run.face == Face::Primary
  {
    let style = face_style(Face::Primary, size, weight, color, anchor);
    area.draw(&Text::new(text, (rounded(x), rounded(y)), style))?;
    return Ok(());
  }

  let total: f64 = runs.iter().map(|run| run.width).sum();
  let left = match anchor.h_pos {
    HPos::Left => x,
    HPos::Center => x - total / 2.0,
    HPos::Right => x - total,
  };
  runs.iter().try_for_each(|run| {
    let style = face_style(run.face, size, weight, color, Pos::new(HPos::Left, anchor.v_pos));
    area.draw(&Text::new(
      run.text.as_str(),
      (rounded(left + run.offset), rounded(y)),
      style,
    ))?;
    Ok(())
  })
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

  fn fallback_bytes(self) -> &'static [u8] {
    match self {
      Self::Regular => FALLBACK_REGULAR,
      Self::Bold => FALLBACK_BOLD,
    }
  }

  fn font_style(self) -> FontStyle {
    match self {
      Self::Regular => FontStyle::Normal,
      Self::Bold => FontStyle::Bold,
    }
  }
}

fn text_runs(text: &str, size: f64, weight: Weight) -> Result<Vec<TextRun>, Report> {
  let size = size
    .to_f32()
    .ok_or_else(|| make_internal_report!("Text size {size} is outside the range of the font scale"))?;
  let faces = ChartFaces::new(weight, size)?;
  let glyphs = text.chars().map(|c| (c, faces.glyph(c))).collect_vec();
  let mut offset = 0.0;
  let runs = glyphs
    .chunk_by(|(_, (left, _)), (_, (right, _))| left == right)
    .map(|run| {
      let face = run.first().map_or(Face::Primary, |&(_, (face, _))| face);
      let font = faces.scaled(face);
      let advances: f32 = run.iter().map(|&(_, (_, glyph))| font.h_advance(glyph)).sum();
      let kerning: f32 = run
        .iter()
        .tuple_windows()
        .map(|(&(_, (_, left)), &(_, (_, right)))| font.kern(left, right))
        .sum();
      let width = f64::from(advances + kerning);
      let text_run = TextRun {
        text: run.iter().map(|&(c, _)| c).collect(),
        face,
        offset,
        width,
      };
      offset += width;
      text_run
    })
    .collect_vec();
  Ok(runs)
}

fn face_style<C: Color>(face: Face, size: f64, weight: Weight, color: &C, anchor: Pos) -> TextStyle<'static> {
  let family = match face {
    Face::Primary => CHART_FONT_FAMILY,
    Face::Fallback => FALLBACK_FONT_FAMILY,
  };
  FontDesc::new(
    FontFamily::Name(family),
    size * PLOTTERS_SIZE_PER_PIXEL,
    weight.font_style(),
  )
  .color(color)
  .pos(anchor)
}

fn rounded(value: f64) -> i32 {
  f64_to_i32_saturating(value.round())
}

fn scaled_font(bytes: &'static [u8], size: f32, weight: Weight) -> Result<PxScaleFont<FontRef<'static>>, Report> {
  let font = FontRef::try_from_slice(bytes)
    .map_err(|_invalid| make_internal_report!("The embedded chart font {weight:?} is not a valid font"))?;
  let units_per_em = font
    .units_per_em()
    .ok_or_else(|| make_internal_report!("The embedded chart font {weight:?} has no units per em"))?;
  let scale = size * font.height_unscaled() / units_per_em;
  Ok(font.into_scaled(scale))
}

struct TextRun {
  text: String,
  face: Face,
  offset: f64,
  width: f64,
}

struct ChartFaces {
  primary: PxScaleFont<FontRef<'static>>,
  fallback: PxScaleFont<FontRef<'static>>,
}

impl ChartFaces {
  fn new(weight: Weight, size: f32) -> Result<Self, Report> {
    Ok(Self {
      primary: scaled_font(weight.font_bytes(), size, weight)?,
      fallback: scaled_font(weight.fallback_bytes(), size, weight)?,
    })
  }

  fn glyph(&self, c: char) -> (Face, GlyphId) {
    let glyph = self.primary.glyph_id(c);
    if glyph.0 == 0 {
      (Face::Fallback, self.fallback.glyph_id(c))
    } else {
      (Face::Primary, glyph)
    }
  }

  fn scaled(&self, face: Face) -> &PxScaleFont<FontRef<'static>> {
    match face {
      Face::Primary => &self.primary,
      Face::Fallback => &self.fallback,
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Face {
  Primary,
  Fallback,
}
