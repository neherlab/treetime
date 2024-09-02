use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::rtt::ClockRegressionResult;
use crate::o;
use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::UTF8_FULL;
use comfy_table::{ContentArrangement, Table};
use crossterm::terminal;
use eyre::Report;
use itertools::{chain, Itertools};
use num_traits::clamp;
use plotters::coord::Shift;
use plotters::prelude::*;
use rgb::RGB8;
use std::path::Path;
use textplots::{Chart, ColorPlot, Shape};

#[cfg(feature = "png")]
use crate::io::file::create_file_or_stdout;
use crate::utils::float_fmt::float_to_significant_digits;
#[cfg(feature = "png")]
use image::{codecs::png::PngEncoder, ColorType, DynamicImage, ImageBuffer, ImageEncoder, Rgb};
use plotters::coord::types::RangedCoordf32;

pub fn write_clock_regression_chart_svg(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: impl AsRef<Path>,
) -> Result<(), Report> {
  let svg = SVGBackend::new(filepath.as_ref(), (1200, 800)).into_drawing_area();
  draw_chart(results, clock_model, &svg)?;
  svg.present()?;
  Ok(())
}

#[cfg(feature = "png")]
pub fn write_clock_regression_chart_png(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: impl AsRef<Path>,
) -> Result<(), Report> {
  let img = write_clock_regression_chart_bitmap(results, clock_model)?;
  let mut f = &mut create_file_or_stdout(filepath)?;
  let encoder = PngEncoder::new(&mut f);
  encoder.write_image(img.as_bytes(), img.width(), img.height(), ColorType::Rgb8)?;
  Ok(())
}

#[cfg(not(feature = "png"))]
pub fn write_clock_regression_chart_png(
  _results: &[ClockRegressionResult],
  _clock_model: &ClockModel,
  _filepath: impl AsRef<Path>,
) -> Result<(), Report> {
  Ok(())
}

#[cfg(feature = "png")]
pub fn write_clock_regression_chart_bitmap(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
) -> Result<DynamicImage, Report> {
  let (width, height) = (1200, 800);
  let mut img: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::new(width, height);
  {
    let bitmap = BitMapBackend::with_buffer(img.as_flat_samples_mut().samples, (width, height)).into_drawing_area();
    draw_chart(results, clock_model, &bitmap)?;
    bitmap.present()?;
  }
  Ok(DynamicImage::ImageRgb8(img))
}

fn draw_chart<'a, DB>(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  drawing_area: &'a DrawingArea<DB, Shift>,
) -> Result<(), Report>
where
  DB: DrawingBackend + 'a,
  <DB as DrawingBackend>::ErrorType: 'static,
{
  let PointsResult {
    norm_points,
    outlier_points,
    line,
    x_min,
    x_max,
    y_min,
    y_max,
  } = gather_points(results, clock_model)?;

  drawing_area.fill(&WHITE)?;
  drawing_area.margin(10, 10, 10, 10);

  let mut chart = ChartBuilder::on(drawing_area)
    .caption(
      "Clock regression",
      ("sans-serif", 40).into_font().style(FontStyle::Bold),
    )
    .x_label_area_size(60)
    .y_label_area_size(60)
    .margin(10)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

  chart
    .configure_mesh()
    .light_line_style(WHITE)
    .bold_line_style(RGBAColor(128, 128, 128, 0.5))
    .x_desc("Date")
    .x_label_formatter(&|&date| float_to_significant_digits(date, 6))
    .y_desc("Div")
    .y_label_formatter(&|&div| float_to_significant_digits(div, 3))
    .axis_desc_style(("sans-serif", 20).into_font().style(FontStyle::Bold))
    .x_label_style(("sans-serif", 16))
    .y_label_style(("sans-serif", 16))
    .draw()?;

  let norm_point_color = RGBColor(8, 232, 140);
  chart
    .draw_series(PointSeries::of_element(
      norm_points,
      4,
      &norm_point_color,
      &|c, s, st| {
        EmptyElement::at(c)
          + Circle::new((0, 0), s, st.filled())
          + Circle::new((0, 0), s, RGBAColor(128, 128, 128, 0.25).stroke_width(1))
      },
    ))?
    .label("Samples (norm)")
    .legend(move |(x, y)| {
      Circle::new(
        (x + 10, y),
        4,
        ShapeStyle {
          color: norm_point_color.to_rgba(),
          filled: true,
          stroke_width: 0,
        },
      )
    });

  let outlier_point_color = RGBColor(255, 105, 97);
  chart
    .draw_series(PointSeries::of_element(
      outlier_points,
      4,
      &outlier_point_color,
      &|c, s, st| {
        EmptyElement::at(c)
          + Circle::new((0, 0), s, st.filled())
          + Circle::new((0, 0), s, RGBAColor(128, 128, 128, 0.25).stroke_width(1))
      },
    ))?
    .label("Samples (outliers)")
    .legend(move |(x, y)| {
      Circle::new(
        (x + 10, y),
        4,
        ShapeStyle {
          color: outlier_point_color.to_rgba(),
          filled: true,
          stroke_width: 0,
        },
      )
    });

  let line_color = RGBColor(8, 140, 232);
  chart
    .draw_series(LineSeries::new(
      line,
      ShapeStyle {
        color: line_color.to_rgba(),
        filled: true,
        stroke_width: 3,
      },
    ))?
    .label(format!("Clock regression: {}", clock_model.equation_str()))
    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], line_color));

  add_legend(
    &mut chart,
    format!("tMRCA: {:}", float_to_significant_digits(clock_model.t_mrca(), 5)),
  )?;
  add_legend(
    &mut chart,
    format!("Rate: {:}", float_to_significant_digits(clock_model.clock_rate(), 3)),
  )?;
  add_legend(
    &mut chart,
    format!(
      "Intercept: {:}",
      float_to_significant_digits(clock_model.intercept(), 4)
    ),
  )?;
  add_legend(
    &mut chart,
    format!("R: {:}", float_to_significant_digits(clock_model.r_val(), 4)),
  )?;
  add_legend(
    &mut chart,
    format!("R²: {:}", float_to_significant_digits(clock_model.r_val().powf(2.0), 4)),
  )?;
  add_legend(&mut chart, format!("χ²: {:.3e}", clock_model.chisq()))?;

  chart
    .configure_series_labels()
    .label_font(("sans-serif", 20).into_font())
    .margin(5)
    .border_style(RGBAColor(128, 128, 128, 0.5))
    .background_style(RGBAColor(200, 200, 200, 0.5))
    .position(SeriesLabelPosition::UpperLeft)
    .draw()?;

  Ok(())
}

fn add_legend<'a, DB>(
  chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf32, RangedCoordf32>>,
  label: impl Into<String>,
) -> Result<(), Report>
where
  DB: DrawingBackend + 'a,
  <DB as DrawingBackend>::ErrorType: 'static,
{
  chart.draw_series(empty_line())?.label(label);
  Ok(())
}

pub fn print_clock_regression_chart(results: &[ClockRegressionResult], clock_model: &ClockModel) -> Result<(), Report> {
  let mut table = Table::new();
  table
    .load_preset(UTF8_FULL)
    .apply_modifier(UTF8_ROUND_CORNERS)
    .apply_modifier(UTF8_SOLID_INNER_BORDERS)
    .set_content_arrangement(ContentArrangement::Dynamic);

  table.add_row([o!("Clock regression"), clock_model.equation_str()]);
  table.add_row([o!("tMRCA"), format!("{:.1}", clock_model.t_mrca())]);
  table.add_row([o!("Rate"), format!("{:.6}", clock_model.clock_rate())]);
  table.add_row([o!("Intercept"), format!("{:.4}", clock_model.intercept())]);
  table.add_row([o!("R"), format!("{:.4}", clock_model.r_val())]);
  table.add_row([o!("R²"), format!("{:.4}", clock_model.r_val().powf(2.0))]);
  table.add_row([o!("χ²"), format!("{:.3e}", clock_model.chisq())]);
  println!("{table}");

  let (width, height) = terminal::size()?;
  let width = clamp(width, 0, 1024) as u32;
  let height = clamp(height, 0, 1024) as u32;

  let PointsResult {
    norm_points,
    outlier_points,
    x_min,
    x_max,
    y_min,
    y_max,
    ..
  } = gather_points(results, clock_model)?;

  let mut chart = Chart::new_with_y_range(width, height, x_min, x_max, y_min, y_max);

  let norm_points = Shape::Points(&norm_points);
  let chart = chart.linecolorplot(&norm_points, RGB8 { r: 8, g: 232, b: 140 });

  let outlier_points = Shape::Points(&outlier_points);
  let chart = chart.linecolorplot(&outlier_points, RGB8 { r: 255, g: 105, b: 97 });

  let line = Box::new(|date: f32| clock_model.div(date as f64) as f32);
  let line = Shape::Continuous(line);
  let chart = chart.linecolorplot(&line, RGB8 { r: 8, g: 140, b: 232 });

  chart.nice();

  Ok(())
}

fn empty_line<DB>() -> LineSeries<DB, (f32, f32)>
where
  DB: DrawingBackend,
{
  LineSeries::<DB, _>::new(Vec::<(f32, f32)>::new(), TRANSPARENT)
}

struct PointsResult {
  norm_points: Vec<(f32, f32)>,
  outlier_points: Vec<(f32, f32)>,
  line: [(f32, f32); 2],
  x_min: f32,
  x_max: f32,
  y_min: f32,
  y_max: f32,
}

fn gather_points(results: &[ClockRegressionResult], clock_model: &ClockModel) -> Result<PointsResult, Report> {
  assert!(!results.is_empty());

  let (norms, outliers): (Vec<_>, Vec<_>) = results.iter().partition(|result| result.is_outlier);

  let norm_points = outliers
    .into_iter()
    .cloned()
    .filter_map(|result| result.date.map(|date| (date as f32, result.div as f32)))
    .collect_vec();

  let outlier_points = norms
    .into_iter()
    .filter_map(|result| result.date.map(|date| (date as f32, result.div as f32)))
    .collect_vec();

  let points = chain!(&norm_points, &outlier_points).copied().collect_vec();

  let (x_min, x_max) = points.iter().map(|(x, _)| *x).minmax().into_option().unwrap();

  let line_y1 = clock_model.div(x_min as f64) as f32;
  let line_y2 = clock_model.div(x_max as f64) as f32;
  let line = [(x_min, line_y1), (x_max, line_y2)];

  let (y_min, y_max) = points
    .iter()
    .map(|(_, y)| *y)
    .chain([line_y1, line_y2])
    .minmax()
    .into_option()
    .unwrap();

  Ok(PointsResult {
    norm_points,
    outlier_points,
    line,
    x_min,
    x_max,
    y_min,
    y_max,
  })
}
