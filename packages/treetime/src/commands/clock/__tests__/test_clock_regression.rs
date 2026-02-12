#[cfg(test)]
mod tests {
  use crate::commands::clock::clock_graph::GraphClock;
  use crate::commands::clock::clock_model::ClockModel;
  use crate::commands::clock::clock_regression::{ClockParams, clock_regression_backward};
  use crate::commands::clock::clock_traits::ClockNode;
  use crate::graph::node::Named;
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::seq::div::{OnlyLeaves, compute_divs};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use std::collections::BTreeMap;

  pub fn compute_naive_rate(dates: &BTreeMap<String, f64>, div: &BTreeMap<String, f64>) -> f64 {
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

    let graph: GraphClock = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let divs = compute_divs(&graph, OnlyLeaves(true));
    let naive_rate = compute_naive_rate(&dates, &divs);

    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().time = Some(dates[&name]);
    }

    clock_regression_backward(&graph, &ClockParams::default());
    let clock = {
      let root = graph.get_exactly_one_root()?;
      let root = root.read_arc().payload().read_arc();
      ClockModel::new(root.clock_set())
    }?;
    assert_ulps_eq!(naive_rate, clock.clock_rate(), max_ulps = 4);

    let options = &ClockParams {
      variance_factor: 1.0,
      variance_offset: 0.0,
      variance_offset_leaf: 1.0,
    };

    clock_regression_backward(&graph, options);
    let clock = {
      let root = graph.get_exactly_one_root()?;
      let root = root.read_arc().payload().read_arc();
      ClockModel::new(root.clock_set())
    }?;
    assert_ulps_eq!(0.007710610998647367, clock.clock_rate(), max_ulps = 4);

    Ok(())
  }
}
