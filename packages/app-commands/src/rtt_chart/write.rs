use crate::rtt_chart::data::RttChartData;
use crate::rtt_chart::frame::RttFrame;
use crate::rtt_chart::layout::{CHART_SIZE, RttLayout, TITLE};
use crate::rtt_chart::render::draw_chart;
use crate::rtt_chart::summary::RttSummary;
use app_output::output_format::ChartFormat;
use eyre::{Report, WrapErr};
use plotters::prelude::{IntoDrawingArea, SVGBackend};
use std::io::Write;
use std::path::Path;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime_charts::convert::{svg_to_pdf, svg_to_png};
use treetime_charts::font::register_chart_fonts;
use treetime_utils::io::file::write_file_with;
use treetime_utils::make_internal_report;

pub fn write_rtt_chart(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: &Path,
) -> Result<(), Report> {
  let format = ChartFormat::from_path(filepath)
    .ok_or_else(|| make_internal_report!("No chart format for '{}' after output planning", filepath.display()))?;
  let bytes = encode_chart(format, rtt_chart_svg(results, clock_model)?)?;
  write_file_with(filepath, |writer| {
    writer
      .write_all(&bytes)
      .wrap_err_with(|| format!("When writing the {format} chart"))
  })
}

pub fn encode_chart(format: ChartFormat, svg: String) -> Result<Vec<u8>, Report> {
  match format {
    ChartFormat::Svg => Ok(svg.into_bytes()),
    ChartFormat::Png => svg_to_png(&svg),
    ChartFormat::Pdf => svg_to_pdf(&svg),
  }
}

pub(crate) fn rtt_chart_svg(results: &[ClockRegressionResult], clock_model: &ClockModel) -> Result<String, Report> {
  register_chart_fonts()?;
  let data = RttChartData::new(results)?;
  let frame = RttFrame::new(&data, clock_model);
  let summary = RttSummary::new(&data, clock_model);
  let layout = RttLayout::new(&data, &frame, &summary)?;
  let mut svg = String::new();
  {
    let area = SVGBackend::with_string(&mut svg, CHART_SIZE).into_drawing_area();
    draw_chart(&area, &data, &frame, &layout)?;
    area.present()?;
  }
  with_accessible_name(&svg, TITLE, &summary.description())
}

fn with_accessible_name(svg: &str, title: &str, description: &str) -> Result<String, Report> {
  let (prefix, root) = svg
    .split_once("<svg")
    .ok_or_else(|| make_internal_report!("The chart library wrote no SVG root element"))?;
  let (attributes, body) = root
    .split_once('>')
    .ok_or_else(|| make_internal_report!("The chart library wrote an unterminated SVG root element"))?;
  Ok(format!(
    "{prefix}<svg role=\"img\"{attributes}>\n<title>{}</title>\n<desc>{}</desc>{body}",
    xml_text(title),
    xml_text(description)
  ))
}

fn xml_text(text: &str) -> String {
  text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}
