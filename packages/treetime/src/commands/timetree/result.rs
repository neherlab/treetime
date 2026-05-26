use crate::clock::clock_model::ClockModel;
use crate::partition::timetree::GraphTimetree;
use crate::timetree::confidence::NodeConfidenceInterval;

#[derive(Debug)]
pub struct TimetreeResult {
  pub graph: GraphTimetree,
  pub clock_model: ClockModel,
  pub confidence_intervals: Option<Vec<NodeConfidenceInterval>>,
}
