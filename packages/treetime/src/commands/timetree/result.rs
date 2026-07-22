use crate::clock::clock_model::ClockModel;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::timetree::GraphTimetree;
use crate::partition::timetree::PartitionTimetreeAllVec;
use crate::timetree::confidence::NodeConfidenceInterval;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_io::dates_csv::DatesMap;

#[derive(Serialize)]
pub struct TimetreeGraphData {
  pub clock_model: ClockModel,
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  pub partitions: PartitionTimetreeAllVec,
  pub dates: Option<DatesMap>,
  pub gtr: Option<GTR>,
  pub model_name: Option<GtrModelName>,
  pub mutation_counts: Option<BTreeMap<treetime_graph::edge::GraphEdgeKey, usize>>,
}

impl TimetreeGraphData {
  pub fn new(
    clock_model: ClockModel,
    confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
    partitions: PartitionTimetreeAllVec,
    dates: Option<DatesMap>,
    gtr: Option<GTR>,
    model_name: Option<GtrModelName>,
    mutation_counts: Option<BTreeMap<treetime_graph::edge::GraphEdgeKey, usize>>,
  ) -> Self {
    Self {
      clock_model,
      confidence_intervals,
      partitions,
      dates,
      gtr,
      model_name,
      mutation_counts,
    }
  }
}

#[derive(Serialize)]
pub struct TimetreeResult {
  #[serde(skip)]
  pub graph: GraphTimetree<TimetreeGraphData>,
}

impl std::ops::Deref for TimetreeResult {
  type Target = TimetreeGraphData;

  fn deref(&self) -> &Self::Target {
    self.graph.data()
  }
}
