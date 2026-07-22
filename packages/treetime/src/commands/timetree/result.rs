use crate::clock::clock_model::ClockModel;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::timetree::GraphTimetree;
use crate::partition::timetree::PartitionTimetreeAllVec;
use crate::timetree::confidence::NodeConfidenceInterval;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_io::dates_csv::DatesMap;

#[allow(clippy::manual_non_exhaustive, clippy::partial_pub_fields)] // The private unit field preserves Graph JSON's `data: null` shape.
#[derive(Serialize)]
#[serde(transparent)]
pub struct TimetreeGraphData {
  marker: (),
  #[serde(skip)]
  pub clock_model: ClockModel,
  #[serde(skip)]
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
  #[serde(skip)]
  pub partitions: PartitionTimetreeAllVec,
  #[serde(skip)]
  pub dates: Option<DatesMap>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub model_name: Option<GtrModelName>,
  #[serde(skip)]
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
      marker: (),
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
