use crate::commands::clock::clock_set::ClockModel;
use crate::commands::clock::rtt::ClockRegressionResult;
use crate::o;
use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::UTF8_FULL;
use comfy_table::{ContentArrangement, Table};
use crossterm::terminal;
use eyre::Report;
use itertools::Itertools;
use num_traits::clamp;
use plotters::coord::Shift;
use plotters::prelude::*;
use rgb::RGB8;
use std::path::Path;
use textplots::{Chart, ColorPlot, Shape};

#[cfg(feature = "png")]
use crate::io::file::create_file_or_stdout;
#[cfg(feature = "png")]
use image::{codecs::png::PngEncoder, ColorType, DynamicImage, ImageBuffer, ImageEncoder, Rgb};

pub fn write_clock_regression_chart_svg(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: impl AsRef<Path>,
) -> Result<(), Report> {
  let svg = SVGBackend::new(filepath.as_ref(), (640, 480)).into_drawing_area();
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
  let (width, height) = (640, 480);
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
    points,
    line,
    x_min,
    x_max,
    y_min,
    y_max,
  } = gather_points(results, clock_model)?;

  drawing_area.fill(&WHITE)?;
  drawing_area.margin(10, 10, 10, 10);

  let mut chart = ChartBuilder::on(drawing_area)
    .caption("Clock regression", ("sans-serif", 20).into_font())
    .x_label_area_size(30)
    .y_label_area_size(60)
    .margin(5)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

  chart
    .configure_mesh()
    .x_desc("Date")
    .y_desc("Div")
    .axis_desc_style(("sans-serif", 12))
    .draw()?;

  let line_color = RGBColor(8, 140, 232);
  chart
    .draw_series(LineSeries::new(line, &line_color))?
    .label(format!("Clock regression: {}", clock_model.equation_str()))
    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], line_color));

  let point_color = RGBColor(255, 105, 97);
  chart
    .draw_series(PointSeries::of_element(points, 2, &point_color, &|c, s, st| {
      EmptyElement::at(c) + Circle::new((0, 0), s, st.filled())
    }))?
    .label("Samples")
    .legend(move |(x, y)| {
      Circle::new(
        (x + 10, y),
        2,
        ShapeStyle {
          color: point_color.to_rgba(),
          filled: true,
          stroke_width: 0,
        },
      )
    });

  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE)
    .position(SeriesLabelPosition::UpperRight)
    .draw()?;

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
  table.add_row([o!("Rate"), format!("{:.4}", clock_model.clock_rate())]);
  table.add_row([o!("Intercept"), format!("{:.4}", clock_model.intercept())]);
  table.add_row([o!("R"), format!("{:.4}", clock_model.r_val())]);
  table.add_row([o!("R²"), format!("{:.4}", clock_model.r_val().powf(2.0))]);
  table.add_row([o!("χ²"), format!("{:.4}", clock_model.chisq())]);
  println!("{table}");

  let (width, height) = terminal::size()?;
  let width = clamp(width, 0, 1024) as u32;
  let height = clamp(height, 0, 1024) as u32;

  let PointsResult {
    points,
    x_min,
    x_max,
    y_min,
    y_max,
    ..
  } = gather_points(results, clock_model)?;

  let points = Shape::Points(&points);

  let line = Box::new(|date: f32| clock_model.div(date as f64) as f32);
  let line = Shape::Continuous(line);

  let mut chart = Chart::new_with_y_range(width, height, x_min, x_max, y_min, y_max);
  let chart = chart.linecolorplot(&line, RGB8 { r: 8, g: 140, b: 232 });
  let chart = chart.linecolorplot(&points, RGB8 { r: 255, g: 105, b: 97 });

  chart.nice();

  Ok(())
}

struct PointsResult {
  points: Vec<(f32, f32)>,
  line: [(f32, f32); 2],
  x_min: f32,
  x_max: f32,
  y_min: f32,
  y_max: f32,
}

fn gather_points(results: &[ClockRegressionResult], clock_model: &ClockModel) -> Result<PointsResult, Report> {
  assert!(!results.is_empty());

  let points = results
    .iter()
    .filter_map(|result| result.date.map(|date| (date as f32, result.div as f32)))
    .collect_vec();

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
    points,
    line,
    x_min,
    x_max,
    y_min,
    y_max,
  })
}
