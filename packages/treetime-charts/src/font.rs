use eyre::Report;
use fontdb::{Database, Source};
use plotters::style::{FontStyle, register_font};
use std::sync::{Arc, LazyLock, OnceLock};
use treetime_utils::make_internal_report;

pub const CHART_FONT_FAMILY: &str = "Lato, sans-serif";
pub const LATO_REGULAR: &[u8] = include_bytes!("../fonts/Lato-Regular.ttf");
pub const LATO_BOLD: &[u8] = include_bytes!("../fonts/Lato-Bold.ttf");

pub(crate) const LATO_FAMILY: &str = "Lato";
const PLOTTERS_FAMILIES: [&str; 2] = [CHART_FONT_FAMILY, "sans-serif"];

static REGISTRATION: OnceLock<Result<(), &'static str>> = OnceLock::new();

pub(crate) static CHART_FONT_DATABASE: LazyLock<Arc<Database>> = LazyLock::new(chart_font_database);

pub fn register_chart_fonts() -> Result<(), Report> {
  let registration = *REGISTRATION.get_or_init(register_lato);
  registration.map_err(|font| make_internal_report!("The chart library rejected the embedded font {font}"))
}

fn register_lato() -> Result<(), &'static str> {
  PLOTTERS_FAMILIES.iter().try_for_each(|family| {
    register_font(family, FontStyle::Normal, LATO_REGULAR).map_err(|_invalid| "Lato Regular")?;
    register_font(family, FontStyle::Bold, LATO_BOLD).map_err(|_invalid| "Lato Bold")
  })
}

fn chart_font_database() -> Arc<Database> {
  let mut database = Database::new();
  database.load_font_source(Source::Binary(Arc::new(LATO_REGULAR)));
  database.load_font_source(Source::Binary(Arc::new(LATO_BOLD)));
  database.set_sans_serif_family(LATO_FAMILY);
  Arc::new(database)
}
