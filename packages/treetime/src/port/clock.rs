use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::{invert_edge, GraphEdge, GraphEdgeKey, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey, Named};
use crate::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use crate::io::nwk::{format_weight, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};
use crate::o;
use crate::port::clock_set::{ClockModel, ClockSet};
use crate::utils::container::get_exactly_one;
use approx::ulps_eq;
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::sync::Arc;

pub type ClockGraph = Graph<ClockNodePayload, ClockEdgePayload, ()>;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct ClockNodePayload {
  pub name: Option<String>,
  pub date: Option<f64>,
  pub div: f64,
  pub is_outlier: bool,
  pub total: ClockSet,
  pub to_parent: ClockSet,
  pub to_children: BTreeMap<String, ClockSet>,
  pub from_children: BTreeMap<String, ClockSet>,
}

impl GraphNode for ClockNodePayload {}

impl NodeFromNwk for ClockNodePayload {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..ClockNodePayload::default()
    })
  }
}

impl NodeToNwk for ClockNodePayload {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mutations: String = "".to_owned(); // TODO: fill mutations
    BTreeMap::from([(o!("mutations"), mutations)])
  }
}

impl Named for ClockNodePayload {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl NodeToGraphviz for ClockNodePayload {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct ClockEdgePayload {
  branch_length: Option<f64>,
}

impl GraphEdge for ClockEdgePayload {}

impl Weighted for ClockEdgePayload {
  fn weight(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_weight(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for ClockEdgePayload {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self { branch_length })
  }
}

impl EdgeToNwk for ClockEdgePayload {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight()
  }
}

impl EdgeToGraphViz for ClockEdgePayload {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight()
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClockOptions {
  pub variance_factor: f64,
  pub variance_offset: f64,
}

fn clock_regression_backward(graph: &ClockGraph, options: &ClockOptions) {
  graph.par_iter_breadth_first_backward(|mut n| {
    n.payload.from_children.clear();
    let q_dest = if n.is_leaf {
      ClockSet::leaf_contribution(n.payload.date)
    } else {
      let mut q_dest = ClockSet::default();
      for (c, e) in n.children {
        let c = c.read_arc();
        let e = e.read_arc();

        let name = c
          .name()
          .expect("Encountered a child node without a name")
          .as_ref()
          .to_owned();

        let edge_len = e.weight().expect("Encountered an edge without a weight");

        let branch_variance = options.variance_factor * edge_len + options.variance_offset;
        let q_from_child = c.to_parent.propagate_averages(edge_len, branch_variance);
        q_dest += &q_from_child;
        n.payload.from_children.insert(name, q_from_child);
      }
      q_dest
    };
    n.payload.to_parent = q_dest;
    GraphTraversalContinuation::Continue
  });
}

fn clock_regression_forward(graph: &ClockGraph, options: &ClockOptions) {
  graph.par_iter_breadth_first_forward(|mut n| {
    let q = &n.payload.to_parent;
    let mut q_tot = q.clone();

    n.payload.div = 0.0;
    n.payload.is_outlier = false; // TODO: calculate this

    if !n.is_root {
      let (p, e) = &n.get_exactly_one_parent().unwrap();
      let p = p.read_arc();
      let e = e.read_arc();

      let name = n
        .payload
        .name()
        .expect("Encountered a node without a name")
        .as_ref()
        .to_owned();

      let branch_length = e.weight().expect("Encountered an edge without a weight");
      let branch_variance = options.variance_factor * branch_length + options.variance_offset;
      q_tot += p.to_children[&name].propagate_averages(branch_length, branch_variance);

      n.payload.div = n.payload.div + branch_length;
    }

    n.payload.to_children = n
      .payload
      .from_children
      .iter()
      .map(|(c, q)| (c.clone(), &q_tot - q))
      .collect();

    n.payload.total = q_tot;

    GraphTraversalContinuation::Continue
  });
}

/// Calculate tip-to-root regression
pub fn run_clock_regression(graph: &ClockGraph, options: &ClockOptions) -> Result<ClockModel, Report> {
  clock_regression_backward(graph, options);
  clock_regression_forward(graph, options);
  let root = graph.get_exactly_one_root()?;
  let clock = root.read_arc().payload().read_arc().total.clock_model()?;
  Ok(clock)
}

pub fn reroot_in_place(graph: &mut ClockGraph, options: &ClockOptions) -> Result<GraphNodeKey, Report> {
  let FindRootResult { edge, split, .. } = find_best_root(graph, options)?;

  let edge_key = edge.expect("Edge is empty when rerooting");
  let edge = graph.get_edge(edge_key).expect("Edge not found");

  let new_root_key = if ulps_eq!(split, 0.0, max_ulps = 5) {
    edge.read_arc().target()
  } else if ulps_eq!(split, 1.0, max_ulps = 5) {
    edge.read_arc().source()
  } else {
    create_new_root_node(graph, edge_key, split)?
  };

  let old_root_key = { graph.get_exactly_one_root()?.read_arc().key() };
  if new_root_key != old_root_key {
    apply_reroot(graph, old_root_key, new_root_key)?;
  }

  // TODO: remove old root node if it's trivial (i.e. having exactly 1 child and 1 parent) and merge dangling edges

  Ok(new_root_key)
}

/// Create new root node by splitting the edge into two
fn create_new_root_node(graph: &mut ClockGraph, edge_key: GraphEdgeKey, split: f64) -> Result<GraphNodeKey, Report> {
  let new_root_key = graph.add_node(ClockNodePayload {
    name: Some("new_root".to_owned()),
    date: None,
    div: 0.0,
    is_outlier: false,
    total: ClockSet::default(),
    to_parent: ClockSet::default(),
    to_children: btreemap! {},
    from_children: btreemap! {},
  });

  let edge = graph.get_edge(edge_key).expect("Edge not found");
  let source_key = edge.read_arc().source();
  let target_key = edge.read_arc().target();
  let branch_length = edge.read_arc().payload().read_arc().weight().unwrap_or_default();

  graph.add_edge(
    source_key,
    new_root_key,
    ClockEdgePayload {
      branch_length: Some(split * branch_length),
    },
  )?;

  graph.add_edge(
    new_root_key,
    target_key,
    ClockEdgePayload {
      branch_length: Some((1.0 - split) * branch_length),
    },
  )?;

  graph.remove_edge(edge_key)?;

  Ok(new_root_key)
}

/// Modify graph topology to make the newly identified root the actual root.
fn apply_reroot(graph: &mut ClockGraph, old_root_key: GraphNodeKey, new_root_key: GraphNodeKey) -> Result<(), Report> {
  // Find paths from the old root to the new desired root
  let paths = graph.path_from_node_to_node(new_root_key, old_root_key)?;

  // Invert every edge on the path from old to new root.
  // This will make the desired new root into an actual root. The old root might no longer be a root.
  for (_, edge) in &paths {
    if let Some(edge) = edge {
      invert_edge(graph, edge);
    }
  }

  // Some bookkeeping
  graph.build()?;
  Ok(())
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FindRootResult {
  edge: Option<GraphEdgeKey>,

  /// The best root might be somewhere half-way on an existing edge.
  /// The position of this best root is the split. If this is 0 or 1, then we reroot on either the parent (source) or
  /// the child (target) of the edge. If this is 0 < x < 1, we put a new node at that point, and reroot on that new node.
  split: f64,

  clock: ClockModel,
}

/// Find the best new root node
fn find_best_root(graph: &ClockGraph, options: &ClockOptions) -> Result<FindRootResult, Report> {
  // Loop over all nodes, pick the one with the lowest chisq, then optimize position along surrounding branches.
  clock_regression_backward(graph, options);
  clock_regression_forward(graph, options);

  let root = graph.get_exactly_one_root()?;

  let clock = root.read_arc().payload().read_arc().total.clock_model()?;
  let mut best_chisq = clock.chisq();
  let mut best_root_node = Arc::clone(&root);
  let mut best_res = FindRootResult {
    edge: None,
    split: 0.0,
    clock,
  };

  // find best node
  for n in graph.get_nodes() {
    let tmp_chisq = n.read_arc().payload().read_arc().total.clock_model()?.chisq();
    if tmp_chisq < best_chisq {
      best_chisq = tmp_chisq;
      best_root_node = Arc::clone(&n);
    }
  }

  let best_root_node = best_root_node.read_arc();

  // Check if someone on parent branch is better
  if !best_root_node.is_root() {
    // TODO: should we loop over parents as well? (just like over children in the loop just below)
    let inbound = best_root_node.inbound();
    let edge = get_exactly_one(inbound).expect("Not implemented: multiple parent nodes");
    let res = find_best_split(graph, *edge, options)?;
    if res.clock.chisq() < best_chisq {
      best_chisq = res.clock.chisq();
      best_res = res;
    }
  }

  // Check if someone on a child branch is better
  for e in best_root_node.outbound() {
    let res = find_best_split(graph, *e, options)?;
    if res.clock.chisq() < best_chisq {
      best_chisq = res.clock.chisq();
      best_res = res;
    }
  }

  Ok(best_res)
}

fn find_best_split(graph: &ClockGraph, edge: GraphEdgeKey, options: &ClockOptions) -> Result<FindRootResult, Report> {
  let edge = graph.get_edge(edge).expect("Edge not found");

  // Get clock data from both ends of the edge.
  let n = graph.get_node(edge.read_arc().target()).expect("Node not found");
  let n_clock = &n.read_arc().payload().read_arc().to_parent;
  let n = n.read_arc().payload().read_arc();
  let n_name = n.name().expect("Encountered node without a name").as_ref().to_owned();

  let p = graph.get_node(edge.read_arc().source()).expect("Node not found");
  let p_clock = &p.read_arc().payload().read_arc().to_children[&n_name];

  // Precalculate branch values for the edge.
  let branch_length = edge
    .read_arc()
    .payload()
    .read_arc()
    .weight()
    .expect("Encountered an edge without a weight");

  let branch_variance = options.variance_factor * branch_length + options.variance_offset;

  // Interrogate different positions along the branch
  let mut best_chisq = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_clock = ClockModel::default();

  // TODO: arbitrary choice for now, should optimize
  for x in Array1::linspace(0.0, 1.0, 11) {
    let Q = p_clock.propagate_averages(branch_length * x, branch_variance * x)
      + n_clock.propagate_averages(branch_length * (1.0 - x), branch_variance * (1.0 - x));
    let clock_model = Q.clock_model()?;
    if clock_model.chisq() < best_chisq {
      best_chisq = clock_model.chisq();
      best_split = x;
      best_clock = clock_model;
    }
  }

  Ok(FindRootResult {
    edge: Some(edge.read_arc().key()),
    split: best_split,
    clock: best_clock,
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::port::div::{calculate_divs, OnlyLeaves};
  use approx::assert_ulps_eq;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;

  pub fn calculate_naive_rate(dates: &BTreeMap<String, f64>, div: &BTreeMap<String, f64>) -> f64 {
    let t: f64 = dates.values().sum();
    let tsq: f64 = dates.values().map(|&x| x * x).sum();
    let dt: f64 = div.iter().map(|(c, div)| div * dates[c]).sum();
    let d: f64 = div.values().sum();
    (dt * 4.0 - d * t) / (tsq * 4.0 - (t * t))
  }

  #[test]
  fn test_calculate_divs_only_leaves() -> Result<(), Report> {
    let graph: ClockGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let actual = calculate_divs(&graph, OnlyLeaves(true));
    let expected = btreemap! {
      o!("A") => 0.20000000298023224,
      o!("B") => 0.30000000447034836,
      o!("C") => 0.2500000037252903,
      o!("D") => 0.16999999806284904,
    };
    assert_eq!(&expected, &actual);
    Ok(())
  }

  #[test]
  fn test_clock_naive_rate() -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let graph: ClockGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let divs = calculate_divs(&graph, OnlyLeaves(true));
    let naive_rate = calculate_naive_rate(&dates, &divs);

    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().date = Some(dates[&name]);
    }

    let options = &mut ClockOptions {
      variance_factor: 0.0,
      variance_offset: 0.0,
    };

    let clock = run_clock_regression(&graph, options)?;
    assert_ulps_eq!(naive_rate, clock.clock_rate(), epsilon = 1e-9);

    options.variance_factor = 1.0;
    let res = find_best_root(&graph, options)?;
    assert_ulps_eq!(0.008095476518345305, res.clock.clock_rate(), epsilon = 1e-9);

    Ok(())
  }
}
