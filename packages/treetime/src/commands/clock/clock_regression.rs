use crate::commands::clock::clock_graph::ClockGraph;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_set::ClockSet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use crate::graph::node::Named;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClockOptions {
  pub variance_factor: f64,
  pub variance_offset: f64,
}

pub fn clock_regression_backward(graph: &ClockGraph, options: &ClockOptions) {
  graph.par_iter_breadth_first_backward(|mut n| {
    n.payload.from_children.clear();
    let q_dest = if n.is_leaf {
      if n.payload.is_outlier {
        ClockSet::outlier_contribution()
      } else {
        ClockSet::leaf_contribution(n.payload.date)
      }
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

pub fn clock_regression_forward(graph: &ClockGraph, options: &ClockOptions) {
  graph.par_iter_breadth_first_forward(|mut n| {
    let q = &n.payload.to_parent;
    let mut q_tot = q.clone();

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
  let root = root.read_arc().payload().read_arc();
  let clock = ClockModel::new(&root.total)?;
  Ok(clock)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::clock::find_best_root::find_best_root;
  use crate::graph::node::Named;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::seq::div::{calculate_divs, OnlyLeaves};
  use approx::assert_ulps_eq;
  use maplit::btreemap;
  use std::collections::BTreeMap;

  pub fn calculate_naive_rate(dates: &BTreeMap<String, f64>, div: &BTreeMap<String, f64>) -> f64 {
    let t: f64 = dates.values().sum();
    let tsq: f64 = dates.values().map(|&x| x * x).sum();
    let dt: f64 = div.iter().map(|(c, div)| div * dates[c]).sum();
    let d: f64 = div.values().sum();
    (dt * 4.0 - d * t) / (tsq * 4.0 - (t * t))
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
    let clock_rate = ClockModel::new(&res.total)?.clock_rate();
    assert_ulps_eq!(0.008095476518345305, clock_rate, epsilon = 1e-9);

    Ok(())
  }
}
