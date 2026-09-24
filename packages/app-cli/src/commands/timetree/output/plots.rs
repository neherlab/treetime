use eyre::Report;
use std::path::Path;
use treetime_graph::graph::Graph;

pub(crate) fn plot_root_to_tip(_graph: &Graph, _out_base: &Path) -> Result<(), Report> {
  todo!("Extract leaf dates and distances, plot with regression line, annotate with R² and outliers")
}

pub(crate) fn plot_time_tree(_graph: &Graph, _out_base: &Path) -> Result<(), Report> {
  todo!("Generate time-scaled tree visualization")
}
