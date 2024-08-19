use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::{GraphEdge, GraphEdgeKey, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::nwk::{EdgeFromNwk, NodeFromNwk};
use crate::port::clock_set::{ClockModel, ClockSet};
use crate::utils::container::get_exactly_one;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::sync::Arc;

pub type ClockGraph = Graph<ClockNodePayload, ClockEdgePayload, ()>;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct ClockNodePayload {
  name: Option<String>,
  date: Option<f64>,
  total: ClockSet,
  to_parent: ClockSet,
  to_children: BTreeMap<String, ClockSet>,
  from_children: BTreeMap<String, ClockSet>,
}

impl NodeFromNwk for ClockNodePayload {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..ClockNodePayload::default()
    })
  }
}

impl GraphNode for ClockNodePayload {}

impl Named for ClockNodePayload {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
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
  fn from_nwk(_: Option<f64>) -> Result<Self, Report> {
    Ok(Self::default())
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClockOptions {
  variance_factor: f64,
  variance_offset: f64,
}

fn clock_regression_backward(graph: &ClockGraph, options: &ClockOptions) {
  graph.par_iter_breadth_first_backward(|mut n| {
    n.payload.from_children.clear();
    let q_dest = if n.is_leaf {
      ClockSet::leaf_contribution(n.payload.date.expect("Encountered a leaf node without a date"))
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

    if !n.is_root {
      let (p, e) = &n.get_exactly_one_parent().unwrap();
      let p = p.read_arc();
      let e = e.read_arc();

      let name = p
        .name()
        .expect("Encountered a parent node without a name")
        .as_ref()
        .to_owned();

      let branch_length = e.weight().expect("Encountered an edge without a weight");
      let branch_variance = options.variance_factor * branch_length + options.variance_offset;
      q_tot += p.to_children[&name].propagate_averages(branch_length, branch_variance);
    }

    // n.payload.to_children.clear();
    // for (c, q) in &n.payload.from_children {
    //   let mut q_to_children = q_tot.clone();
    //   q_to_children -= q;
    //   n.payload.to_children.insert(c.clone(), q_to_children);
    // }

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

#[derive(Debug, Serialize, Deserialize)]
struct FindRootResult {
  edge: Option<GraphEdgeKey>,
  split: f64,
  clock: ClockModel,
}

pub fn find_best_root(graph: &ClockGraph, options: &ClockOptions) -> Result<FindRootResult, Report> {
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
      best_root_node = Arc::clone(n);
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
  use approx::assert_ulps_eq;
  use maplit::btreemap;

  #[test]
  fn test_clock_naive_rate() -> Result<(), Report> {
    let dates = btreemap! {
        o!("A") => 2013.0,
        o!("B") => 2022.0,
        o!("C") => 2017.0,
        o!("D") => 2005.0,
    };

    let div = btreemap! {
        o!("A") => 0.2,
        o!("B") => 0.3,
        o!("C") => 0.25,
        o!("D") => 0.17,
    };

    let t: f64 = dates.values().sum();
    let tsq: f64 = dates.values().map(|&x| x * x).sum();
    let dt: f64 = div.iter().map(|(c, div)| div * dates[c]).sum();
    let d: f64 = div.values().sum();
    let naive_rate = (dt * 4.0 - d * t) / (tsq * 4.0 - (t * t));

    let graph: ClockGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().date = Some(dates[&name]);
    }

    let options = &mut ClockOptions {
      variance_factor: 0.0,
      variance_offset: 0.0,
    };

    clock_regression_backward(&graph, options);
    clock_regression_forward(&graph, options);

    let root = graph.get_exactly_one_root()?;
    let clock = root.read_arc().payload().read_arc().total.clock_model()?;

    assert_ulps_eq!(naive_rate, clock.rate());

    options.variance_factor = 1.0;
    let res = find_best_root(&graph, options)?;
    assert_ulps_eq!(0.008095476518345305, res.clock.rate());

    Ok(())
  }
}
