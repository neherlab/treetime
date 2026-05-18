use crate::partition::traits::PartitionTimetreeAll;
use crate::partition::payload::timetree::EdgeTimetree;
use crate::partition::payload::timetree::NodeTimetree;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::graph::Graph;

pub type GraphTimetree = Graph<NodeTimetree, EdgeTimetree, ()>;
pub type PartitionTimetreeAllVec = Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>>;
