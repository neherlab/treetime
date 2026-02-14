use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::graph::Graph;

pub type GraphTimetree = Graph<NodeTimetree, EdgeTimetree, ()>;
pub type PartitionTimetreeAllVec = Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>>;
