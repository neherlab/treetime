use crate::font::{CHART_FONT_DATABASE, LATO_FAMILY};
use eyre::{Report, WrapErr};
use krilla::Document;
use krilla::geom::Size;
use krilla::page::PageSettings;
use krilla_svg::{SurfaceExt, SvgSettings};
use std::sync::Arc;
use tiny_skia::{Pixmap, Transform};
use treetime_utils::{make_internal_report, make_report};
use usvg::{Options, Tree};

const PDF_POINTS_PER_PIXEL: f32 = 0.75;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChartFormat {
  Svg,
  Png,
  Pdf,
}

impl ChartFormat {
  pub fn encode(self, svg: &str) -> Result<Vec<u8>, Report> {
    match self {
      Self::Svg => Ok(svg.as_bytes().to_vec()),
      Self::Png => svg_to_png(svg),
      Self::Pdf => svg_to_pdf(svg),
    }
  }
}

pub fn svg_to_png(svg: &str) -> Result<Vec<u8>, Report> {
  let tree = parse_svg(svg)?;
  let size = tree.size().to_int_size();
  let mut pixmap = Pixmap::new(size.width(), size.height()).ok_or_else(|| {
    make_report!(
      "Cannot create a {}x{} pixel image of the chart",
      size.width(),
      size.height()
    )
  })?;
  resvg::render(&tree, Transform::default(), &mut pixmap.as_mut());
  pixmap.encode_png().wrap_err("When encoding the chart as PNG")
}

pub fn svg_to_pdf(svg: &str) -> Result<Vec<u8>, Report> {
  let tree = parse_svg(svg)?;
  let svg_size = tree.size();
  let page_size = Size::from_wh(
    svg_size.width() * PDF_POINTS_PER_PIXEL,
    svg_size.height() * PDF_POINTS_PER_PIXEL,
  )
  .ok_or_else(|| make_report!("Cannot create a PDF page for a chart of size {svg_size:?}"))?;

  let mut document = Document::new();
  let mut page = document.start_page_with(PageSettings::new(page_size));
  let mut surface = page.surface();
  surface
    .draw_svg(&tree, page_size, SvgSettings::default())
    .ok_or_else(|| make_internal_report!("Cannot draw the chart into the PDF page"))?;
  surface.finish();
  page.finish();
  document.finish().wrap_err("When writing the chart as PDF")
}

fn parse_svg(svg: &str) -> Result<Tree, Report> {
  let options = Options {
    font_family: LATO_FAMILY.to_owned(),
    fontdb: Arc::clone(&CHART_FONT_DATABASE),
    ..Options::default()
  };
  Tree::from_str(svg, &options).wrap_err("When reading the SVG chart")
}
