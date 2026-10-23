use crate::rtt_chart::data::RttChartData;
use crate::rtt_chart::frame::{RootPlacement, RttFrame};
use crate::rtt_chart::summary::{RttSummary, Statistic};
use eyre::Report;
use itertools::Itertools;
use plotters::style::RGBColor;
use std::iter::once;
use treetime_charts::palette::{ACCENT, INK, MUTED};
use treetime_charts::scale::LinearScale;
use treetime_charts::text::{Weight, text_width};
use treetime_charts::ticks::{Tick, date_ticks, linear_ticks};
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

pub(crate) const CHART_SIZE: (u32, u32) = (1200, 800);

pub(crate) const TITLE: &str = "Root-to-tip regression";

const Y_AXIS_TITLE: &str = "Root-to-tip divergence (substitutions per site)";

const X_AXIS_TITLE: &str = "Sampling date";

const MARGIN_X: f64 = 40.0;
const MARGIN_TOP: f64 = 36.0;
const MARGIN_BOTTOM: f64 = 30.0;

const TITLE_SIZE: f64 = 28.0;
const SUBTITLE_SIZE: f64 = 14.0;
const SUBTITLE_TOP: f64 = MARGIN_TOP + TITLE_SIZE + 12.0;
const SUBTITLE_LINE_HEIGHT: f64 = 22.0;

const STAT_LABEL_SIZE: f64 = 12.0;
const STAT_VALUE_SIZE: f64 = 24.0;
const STAT_DETAIL_SIZE: f64 = 12.0;
const STAT_LABEL_TOP: f64 = MARGIN_TOP + 2.0;
const STAT_VALUE_TOP: f64 = STAT_LABEL_TOP + STAT_LABEL_SIZE + 7.0;
const STAT_DETAIL_TOP: f64 = STAT_VALUE_TOP + STAT_VALUE_SIZE + 6.0;
const STAT_COLUMN_GAP: f64 = 44.0;
const SUBTITLE_STATS_GAP: f64 = 32.0;

const ROW_GAP: f64 = 28.0;
const ROW_HEIGHT: f64 = 26.0;
const LEGEND_SIZE: f64 = 14.0;
const LEGEND_SYMBOL_WIDTH: f64 = 24.0;
const LEGEND_SYMBOL_GAP: f64 = 8.0;
const LEGEND_ITEM_GAP: f64 = 28.0;
const LEGEND_TITLE_GAP: f64 = 40.0;

const AXIS_TITLE_SIZE: f64 = 14.0;
const PLOT_TOP_GAP: f64 = 18.0;
const TICK_SIZE: f64 = 12.0;
const TICK_GAP: f64 = 10.0;
const X_TITLE_GAP: f64 = 12.0;
const X_TICK_SPACING: f64 = 110.0;
const Y_TICK_SPACING: f64 = 64.0;
const MIN_TICK_COUNT: usize = 2;

const ROOT_LABEL_SIZE: f64 = 12.0;
const ROOT_LABEL_OFFSET: f64 = 10.0;
const ROOT_LABEL_GAP: f64 = 7.0;

const CAP_CENTER: f64 = 0.4;

const MARKER_RADII: [(usize, i32); 3] = [(150, 5), (600, 4), (3000, 3)];
const SMALLEST_MARKER_RADIUS: i32 = 2;

pub(crate) struct RttLayout {
  pub labels: Vec<Label>,
  pub legend: Vec<LegendEntry>,
  pub plot: Rect,
  pub dates: LinearScale,
  pub divergences: LinearScale,
  pub date_ticks: Vec<f64>,
  pub divergence_ticks: Vec<f64>,
  pub marker_radius: i32,
}

impl RttLayout {
  pub(crate) fn new(data: &RttChartData, frame: &RttFrame, summary: &RttSummary) -> Result<Self, Report> {
    let (width, height) = (f64::from(CHART_SIZE.0), f64::from(CHART_SIZE.1));
    let mut labels = vec![Label::new(
      TITLE,
      (MARGIN_X, MARGIN_TOP),
      TextSpec::new(TITLE_SIZE, Weight::Bold, INK),
      Align::Left,
    )?];

    let stats = statistic_labels(summary, width - MARGIN_X)?;
    let stats_left = stats.iter().map(|label| label.bounds().left).fold(width, f64::min);
    let stats_bottom = stats.iter().map(|label| label.bounds().bottom).fold(0.0, f64::max);
    labels.extend(stats);

    let subtitles = subtitle_labels(summary, SUBTITLE_TOP)?;
    let subtitle_right = subtitles.iter().map(|label| label.bounds().right).fold(0.0, f64::max);
    let subtitles = if subtitle_right + SUBTITLE_STATS_GAP > stats_left {
      subtitle_labels(summary, stats_bottom + SUBTITLE_LINE_HEIGHT - SUBTITLE_SIZE)?
    } else {
      subtitles
    };
    labels.extend(subtitles);
    let header_bottom = labels.iter().map(|label| label.bounds().bottom).fold(0.0, f64::max);

    let legend_items = legend_items(data, frame);
    let legend_widths = legend_items
      .iter()
      .map(|(_, text)| Ok(LEGEND_SYMBOL_WIDTH + LEGEND_SYMBOL_GAP + text_width(text, LEGEND_SIZE, Weight::Regular)?))
      .collect::<Result<Vec<_>, Report>>()?;
    let legend_width =
      legend_widths.iter().sum::<f64>() + LEGEND_ITEM_GAP * legend_widths.len().saturating_sub(1).to_f64_exact();
    let y_title_width = text_width(Y_AXIS_TITLE, AXIS_TITLE_SIZE, Weight::Regular)?;
    let row_top = header_bottom + ROW_GAP;
    let y_title_top = if y_title_width + LEGEND_TITLE_GAP + legend_width <= width - 2.0 * MARGIN_X {
      row_top
    } else {
      row_top + ROW_HEIGHT
    };
    let legend = legend_entries(&legend_items, &legend_widths, width - MARGIN_X - legend_width, row_top)?;
    labels.push(Label::new(
      Y_AXIS_TITLE,
      (MARGIN_X, y_title_top),
      TextSpec::new(AXIS_TITLE_SIZE, Weight::Regular, MUTED),
      Align::Left,
    )?);

    let x_title_top = height - MARGIN_BOTTOM - AXIS_TITLE_SIZE;
    let x_tick_top = x_title_top - X_TITLE_GAP - TICK_SIZE;
    let plot_top = y_title_top + AXIS_TITLE_SIZE + PLOT_TOP_GAP;
    let plot_bottom = x_tick_top - TICK_GAP;

    let divergence_ticks = linear_ticks(
      frame.divergences.0,
      frame.divergences.1,
      tick_count(plot_bottom - plot_top, Y_TICK_SPACING),
    );
    let tick_spec = TextSpec::new(TICK_SIZE, Weight::Regular, MUTED);
    let tick_label_width = divergence_ticks
      .iter()
      .map(|tick| text_width(&tick.label, TICK_SIZE, Weight::Regular))
      .fold_ok(0.0, f64::max)?;
    let plot = Rect {
      left: MARGIN_X + tick_label_width + TICK_GAP,
      top: plot_top,
      right: width - MARGIN_X,
      bottom: plot_bottom,
    };
    let dates = LinearScale::new(frame.dates, (plot.left, plot.right));
    let divergences = LinearScale::new(frame.divergences, (plot.bottom, plot.top));
    let date_ticks = date_ticks(frame.dates.0, frame.dates.1, tick_count(plot.width(), X_TICK_SPACING))?;

    for Tick { value, label } in &divergence_ticks {
      let y = divergences.map(*value) - CAP_CENTER * TICK_SIZE;
      labels.push(Label::new(label, (plot.left - TICK_GAP, y), tick_spec, Align::Right)?);
    }
    for Tick { value, label } in &date_ticks {
      labels.push(Label::new(
        label,
        (dates.map(*value), x_tick_top),
        tick_spec,
        Align::Center,
      )?);
    }
    labels.push(Label::new(
      X_AXIS_TITLE,
      (f64::midpoint(plot.left, plot.right), x_title_top),
      TextSpec::new(AXIS_TITLE_SIZE, Weight::Regular, MUTED),
      Align::Center,
    )?);
    labels.extend(root_label(frame, summary, &plot, &dates, &divergences, x_title_top)?);

    Ok(Self {
      labels,
      legend,
      plot,
      dates,
      divergences,
      date_ticks: date_ticks.iter().map(|tick| tick.value).collect(),
      divergence_ticks: divergence_ticks.iter().map(|tick| tick.value).collect(),
      marker_radius: marker_radius(data.dated().count()),
    })
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Rect {
  pub left: f64,
  pub top: f64,
  pub right: f64,
  pub bottom: f64,
}

impl Rect {
  pub(crate) fn width(&self) -> f64 {
    self.right - self.left
  }

  pub(crate) fn center_y(&self) -> f64 {
    f64::midpoint(self.top, self.bottom)
  }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Label {
  pub text: String,
  pub position: (f64, f64),
  pub spec: TextSpec,
  pub align: Align,
  pub role: LabelRole,
  pub width: f64,
}

impl Label {
  fn new(text: &str, position: (f64, f64), spec: TextSpec, align: Align) -> Result<Self, Report> {
    Ok(Self {
      text: text.to_owned(),
      position,
      spec,
      align,
      role: LabelRole::Text,
      width: text_width(text, spec.size, spec.weight)?,
    })
  }

  fn with_role(self, role: LabelRole) -> Self {
    Self { role, ..self }
  }

  pub(crate) fn bounds(&self) -> Rect {
    let (x, top) = self.position;
    let left = match self.align {
      Align::Left => x,
      Align::Center => x - self.width / 2.0,
      Align::Right => x - self.width,
    };
    Rect {
      left,
      top,
      right: left + self.width,
      bottom: top + self.spec.size,
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Align {
  Left,
  Center,
  Right,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LabelRole {
  Text,
  RootMarker,
  RootNote,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct TextSpec {
  pub size: f64,
  pub weight: Weight,
  pub color: RGBColor,
}

impl TextSpec {
  const fn new(size: f64, weight: Weight, color: RGBColor) -> Self {
    Self { size, weight, color }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct LegendEntry {
  pub symbol: LegendSymbol,
  pub symbol_box: Rect,
  pub label: Label,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LegendSymbol {
  Sample,
  Excluded,
  Fit,
  Extrapolation,
}

fn subtitle_labels(summary: &RttSummary, top: f64) -> Result<Vec<Label>, Report> {
  once(summary.samples.as_str())
    .chain(summary.notes.iter().map(String::as_str))
    .enumerate()
    .map(|(line, text)| {
      Label::new(
        text,
        (MARGIN_X, top + line.to_f64_exact() * SUBTITLE_LINE_HEIGHT),
        TextSpec::new(SUBTITLE_SIZE, Weight::Regular, MUTED),
        Align::Left,
      )
    })
    .collect()
}

fn statistic_labels(summary: &RttSummary, right: f64) -> Result<Vec<Label>, Report> {
  let statistics = summary.statistics().collect_vec();
  let mut column_right = right;
  let mut columns = vec![];
  for statistic in statistics.into_iter().rev() {
    let column = statistic_column(statistic, column_right)?;
    column_right = column
      .iter()
      .map(|label| label.bounds().left)
      .fold(column_right, f64::min)
      - STAT_COLUMN_GAP;
    columns.push(column);
  }
  Ok(columns.into_iter().rev().flatten().collect())
}

fn statistic_column(statistic: &Statistic, right: f64) -> Result<Vec<Label>, Report> {
  let rows = [
    Some((
      statistic.label.as_str(),
      STAT_LABEL_TOP,
      TextSpec::new(STAT_LABEL_SIZE, Weight::Regular, MUTED),
    )),
    Some((
      statistic.value.as_str(),
      STAT_VALUE_TOP,
      TextSpec::new(STAT_VALUE_SIZE, Weight::Bold, INK),
    )),
    statistic.detail.as_deref().map(|detail| {
      (
        detail,
        STAT_DETAIL_TOP,
        TextSpec::new(STAT_DETAIL_SIZE, Weight::Regular, MUTED),
      )
    }),
  ];
  let rows = rows.into_iter().flatten().collect_vec();
  let column_width = rows
    .iter()
    .map(|(text, _, spec)| text_width(text, spec.size, spec.weight))
    .fold_ok(0.0, f64::max)?;
  rows
    .into_iter()
    .map(|(text, top, spec)| Label::new(text, (right - column_width, top), spec, Align::Left))
    .collect()
}

fn legend_items(data: &RttChartData, frame: &RttFrame) -> Vec<(LegendSymbol, &'static str)> {
  [
    (!data.fitted.is_empty()).then_some((LegendSymbol::Sample, "Samples in the fit")),
    (!data.excluded.is_empty()).then_some((LegendSymbol::Excluded, "Excluded by the clock filter")),
    Some((LegendSymbol::Fit, "Least-squares fit")),
    frame
      .extrapolation
      .is_some()
      .then_some((LegendSymbol::Extrapolation, "Extrapolation to the root")),
  ]
  .into_iter()
  .flatten()
  .collect()
}

fn legend_entries(
  items: &[(LegendSymbol, &str)],
  widths: &[f64],
  left: f64,
  top: f64,
) -> Result<Vec<LegendEntry>, Report> {
  let center_y = top + CAP_CENTER * LEGEND_SIZE;
  let lefts = widths.iter().scan(left, |next_left, width| {
    let item_left = *next_left;
    *next_left += width + LEGEND_ITEM_GAP;
    Some(item_left)
  });
  items
    .iter()
    .zip(lefts)
    .map(|(&(symbol, text), symbol_left)| {
      Ok(LegendEntry {
        symbol,
        symbol_box: Rect {
          left: symbol_left,
          top: center_y - LEGEND_SIZE / 2.0,
          right: symbol_left + LEGEND_SYMBOL_WIDTH,
          bottom: center_y + LEGEND_SIZE / 2.0,
        },
        label: Label::new(
          text,
          (symbol_left + LEGEND_SYMBOL_WIDTH + LEGEND_SYMBOL_GAP, top),
          TextSpec::new(LEGEND_SIZE, Weight::Regular, INK),
          Align::Left,
        )?,
      })
    })
    .collect()
}

fn root_label(
  frame: &RttFrame,
  summary: &RttSummary,
  plot: &Rect,
  dates: &LinearScale,
  divergences: &LinearScale,
  x_title_top: f64,
) -> Result<Option<Label>, Report> {
  let value = &summary.root_date.value;
  let spec = TextSpec::new(ROOT_LABEL_SIZE, Weight::Bold, ACCENT);
  match frame.root {
    RootPlacement::OnChart(root) => {
      let (x, y) = (dates.map(root), divergences.map(0.0));
      let text = format!("Root {value}");
      let below = y + ROOT_LABEL_GAP;
      let top = if below + ROOT_LABEL_SIZE <= plot.bottom {
        below
      } else {
        y - ROOT_LABEL_GAP - ROOT_LABEL_SIZE
      };
      let label = Label::new(&text, (x + ROOT_LABEL_OFFSET, top), spec, Align::Left)?;
      let label = if label.bounds().right > plot.right {
        Label::new(&text, (x - ROOT_LABEL_OFFSET, top), spec, Align::Right)?
      } else {
        label
      };
      Ok(Some(label.with_role(LabelRole::RootMarker)))
    },
    RootPlacement::BeforeChart(_) => Ok(Some(
      Label::new(
        &format!("\u{2190} Root date {value}"),
        (plot.left, x_title_top),
        spec,
        Align::Left,
      )?
      .with_role(LabelRole::RootNote),
    )),
    RootPlacement::AfterChart(_) => Ok(Some(
      Label::new(
        &format!("Root date {value} \u{2192}"),
        (plot.right, x_title_top),
        spec,
        Align::Right,
      )?
      .with_role(LabelRole::RootNote),
    )),
    RootPlacement::Undefined => Ok(None),
  }
}

fn tick_count(length: f64, spacing: f64) -> usize {
  f64_to_usize_saturating((length / spacing).floor()).max(MIN_TICK_COUNT)
}

fn marker_radius(count: usize) -> i32 {
  MARKER_RADII
    .iter()
    .find(|&&(max_count, _)| count <= max_count)
    .map_or(SMALLEST_MARKER_RADIUS, |&(_, radius)| radius)
}
