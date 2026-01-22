use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::graph::graph::Graph;
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use parking_lot::RwLock;
use std::sync::Arc;

pub type GraphTimetree = Graph<NodeTimetree, EdgeTimetree, ()>;
pub type PartitionTimetreeAllVec = Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>>;
