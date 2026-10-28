use crate::rtt_chart::data::{RttChartData, RttPoint};
use crate::rtt_chart::frame::{RootPlacement, RttFrame};
use crate::rtt_chart::layout::{Align, Label, LegendEntry, LegendSymbol, RttLayout};
use eyre::Report;
use plotters::coord::Shift;
use plotters::element::{Circle, DashedPathElement, PathElement, Text};
use plotters::prelude::{DrawingArea, DrawingBackend};
use plotters::style::Color;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use treetime_charts::palette::{AXIS, FAULT, FIT, GRID, PAPER, SERIES};
use treetime_charts::text::text_style;
use treetime_utils::cast::f64_to_i32_saturating;

const SAMPLE_FILL_OPACITY: f64 = 0.42;
const SAMPLE_STROKE_OPACITY: f64 = 0.9;
const LINE_CASING_OPACITY: f64 = 0.75;
const LINE_WIDTH: u32 = 2;
const LINE_CASING_WIDTH: u32 = 6;
const DASH_LENGTH: i32 = 8;
const DASH_GAP: i32 = 6;
const ROOT_MARKER_RADIUS: i32 = 5;
const ROOT_MARKER_HALO: i32 = 3;

pub(crate) fn draw_chart<DB>(
  area: &DrawingArea<DB, Shift>,
  data: &RttChartData,
  frame: &RttFrame,
  layout: &RttLayout,
) -> Result<(), Report>
where
  DB: DrawingBackend,
  DB::ErrorType: 'static,
{
  area.fill(&PAPER)?;
  draw_grid(area, frame, layout)?;

  let pixel = |&RttPoint { date, div }: &RttPoint| (layout.dates.pixel(date), layout.divergences.pixel(div));
  for point in &data.fitted {
    draw_sample(area, pixel(point), layout.marker_radius)?;
  }
  for point in &data.excluded {
    draw_excluded(area, pixel(point), layout.marker_radius)?;
  }

  let line_pixel = |&(date, div): &(f64, f64)| (layout.dates.pixel(date), layout.divergences.pixel(div));
  if let Some(extrapolation) = &frame.extrapolation {
    area.draw(&DashedPathElement::new(
      extrapolation.iter().map(line_pixel),
      DASH_LENGTH,
      DASH_GAP,
      FIT.stroke_width(LINE_WIDTH),
    ))?;
  }
  let fit = frame.fit.iter().map(line_pixel).collect::<Vec<_>>();
  area.draw(&PathElement::new(
    fit.clone(),
    PAPER.mix(LINE_CASING_OPACITY).stroke_width(LINE_CASING_WIDTH),
  ))?;
  area.draw(&PathElement::new(fit, FIT.stroke_width(LINE_WIDTH)))?;

  if let RootPlacement::OnChart(root) = frame.root {
    let center = (layout.dates.pixel(root), layout.divergences.pixel(0.0));
    area.draw(&Circle::new(
      center,
      ROOT_MARKER_RADIUS + ROOT_MARKER_HALO,
      PAPER.filled(),
    ))?;
    area.draw(&Circle::new(center, ROOT_MARKER_RADIUS, FIT.filled()))?;
  }

  for entry in &layout.legend {
    draw_legend_entry(area, entry, layout.marker_radius)?;
  }
  for label in &layout.labels {
    draw_label(area, label)?;
  }
  Ok(())
}

fn draw_grid<DB>(area: &DrawingArea<DB, Shift>, frame: &RttFrame, layout: &RttLayout) -> Result<(), Report>
where
  DB: DrawingBackend,
  DB::ErrorType: 'static,
{
  let plot = &layout.plot;
  let (left, right) = (rounded(plot.left), rounded(plot.right));
  let (top, bottom) = (rounded(plot.top), rounded(plot.bottom));
  for &date in &layout.date_ticks {
    let x = layout.dates.pixel(date);
    area.draw(&PathElement::new(vec![(x, top), (x, bottom)], GRID.stroke_width(1)))?;
  }
  for &div in &layout.divergence_ticks {
    let y = layout.divergences.pixel(div);
    area.draw(&PathElement::new(vec![(left, y), (right, y)], GRID.stroke_width(1)))?;
  }
  let (low, high) = frame.divergences;
  if low < 0.0 && high > 0.0 {
    let y = layout.divergences.pixel(0.0);
    area.draw(&PathElement::new(vec![(left, y), (right, y)], AXIS.stroke_width(1)))?;
  }
  area.draw(&PathElement::new(
    vec![(left, bottom), (right, bottom)],
    AXIS.stroke_width(1),
  ))?;
  Ok(())
}

fn draw_sample<DB>(area: &DrawingArea<DB, Shift>, center: (i32, i32), radius: i32) -> Result<(), Report>
where
  DB: DrawingBackend,
  DB::ErrorType: 'static,
{
  area.draw(&Circle::new(center, radius, SERIES.mix(SAMPLE_FILL_OPACITY).filled()))?;
  area.draw(&Circle::new(
    center,
    radius,
    SERIES.mix(SAMPLE_STROKE_OPACITY).stroke_width(1),
  ))?;
  Ok(())
}

fn draw_excluded<DB>(area: &DrawingArea<DB, Shift>, center: (i32, i32), radius: i32) -> Result<(), Report>
where
  DB: DrawingBackend,
  DB::ErrorType: 'static,
{
  area.draw(&Circle::new(center, radius + 1, FAULT.stroke_width(2)))?;
  Ok(())
}

fn draw_legend_entry<DB>(area: &DrawingArea<DB, Shift>, entry: &LegendEntry, radius: i32) -> Result<(), Report>
where
  DB: DrawingBackend,
  DB::ErrorType: 'static,
{
  let symbol = &entry.symbol_box;
  let y = rounded(symbol.center_y());
  let center = (rounded(f64::midpoint(symbol.left, symbol.right)), y);
  let ends = vec![(rounded(symbol.left), y), (rounded(symbol.right), y)];
  match entry.symbol {
    LegendSymbol::Sample => draw_sample(area, center, radius)?,
    LegendSymbol::Excluded => draw_excluded(area, center, radius)?,
    LegendSymbol::Fit => area.draw(&PathElement::new(ends, FIT.stroke_width(LINE_WIDTH)))?,
    LegendSymbol::Extrapolation => area.draw(&DashedPathElement::new(
      ends,
      DASH_LENGTH,
      DASH_GAP,
      FIT.stroke_width(LINE_WIDTH),
    ))?,
  }
  draw_label(area, &entry.label)
}

fn draw_label<DB>(area: &DrawingArea<DB, Shift>, label: &Label) -> Result<(), Report>
where
  DB: DrawingBackend,
  DB::ErrorType: 'static,
{
  let align = match label.align {
    Align::Left => HPos::Left,
    Align::Center => HPos::Center,
    Align::Right => HPos::Right,
  };
  let style = text_style(
    label.spec.size,
    label.spec.weight,
    &label.spec.color,
    Pos::new(align, VPos::Top),
  );
  let (x, y) = label.position;
  area.draw(&Text::new(label.text.as_str(), (rounded(x), rounded(y)), style))?;
  Ok(())
}

fn rounded(value: f64) -> i32 {
  f64_to_i32_saturating(value.round())
}
