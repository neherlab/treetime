use eyre::Report;
use fontdb::{Database, Source};
use plotters::style::{FontStyle, register_font};
use std::sync::{Arc, LazyLock, OnceLock};
use treetime_utils::make_internal_report;

pub const CHART_FONT_FAMILY: &str = "IBM Plex Sans, sans-serif";
pub const PLEX_SANS_REGULAR: &[u8] = include_bytes!("../fonts/IBMPlexSans-Regular.ttf");
pub const PLEX_SANS_BOLD: &[u8] = include_bytes!("../fonts/IBMPlexSans-Bold.ttf");

pub(crate) const PLEX_SANS_FAMILY: &str = "IBM Plex Sans";
const PLOTTERS_FAMILIES: [&str; 2] = [CHART_FONT_FAMILY, "sans-serif"];

static REGISTRATION: OnceLock<Result<(), &'static str>> = OnceLock::new();

pub(crate) static CHART_FONT_DATABASE: LazyLock<Arc<Database>> = LazyLock::new(chart_font_database);

pub fn register_chart_fonts() -> Result<(), Report> {
  let registration = *REGISTRATION.get_or_init(register_plex_sans);
  registration.map_err(|font| make_internal_report!("The chart library rejected the embedded font {font}"))
}

fn register_plex_sans() -> Result<(), &'static str> {
  PLOTTERS_FAMILIES.iter().try_for_each(|family| {
    register_font(family, FontStyle::Normal, PLEX_SANS_REGULAR).map_err(|_invalid| "IBM Plex Sans Regular")?;
    register_font(family, FontStyle::Bold, PLEX_SANS_BOLD).map_err(|_invalid| "IBM Plex Sans Bold")
  })
}

fn chart_font_database() -> Arc<Database> {
  let mut database = Database::new();
  database.load_font_source(Source::Binary(Arc::new(PLEX_SANS_REGULAR)));
  database.load_font_source(Source::Binary(Arc::new(PLEX_SANS_BOLD)));
  database.set_sans_serif_family(PLEX_SANS_FAMILY);
  Arc::new(database)
}
