use crate::commands::timetree::timetree_args::TreetimeTimetreeArgs;
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum DateConstraint {
  Point(f64),
  Range { start: f64, end: f64 }, // TODO: enforce start <= end when constructing range constraints.
  Distribution(Arc<Distribution>),
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DateConstraintSet {
  pub per_node: BTreeMap<GraphNodeKey, DateConstraint>,
}

pub fn load_date_constraints<N, E, D>(
  args: &TreetimeTimetreeArgs,
  graph: &Graph<N, E, D>,
) -> Result<DateConstraintSet, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let _ = (args, graph);
  // TODO: parse dates/metadata files and attach constraints to matching nodes.
  todo!()
}

pub fn merge_constraint_sets(
  base: DateConstraintSet,
  supplemental: DateConstraintSet,
) -> Result<DateConstraintSet, Report> {
  drop((base, supplemental));
  // TODO: define precedence rules when multiple sources specify the same node's constraint.
  todo!()
}

pub fn constraint_for(constraints: &DateConstraintSet, key: GraphNodeKey) -> Option<&DateConstraint> {
  constraints.per_node.get(&key)
}
