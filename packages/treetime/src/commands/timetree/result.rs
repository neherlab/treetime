use crate::clock::clock_model::ClockModel;
use crate::partition::timetree::GraphTimetree;
use crate::timetree::confidence::NodeConfidenceInterval;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct TimetreeResult {
  #[serde(skip)]
  pub graph: GraphTimetree,
  pub clock_model: ClockModel,
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
}
