#[cfg(test)]
mod tests {
  use crate::commands::timetree::inference::runner::create_branch_distributions_input_mode;
  use crate::commands::timetree::timetree_traits::TimetreeEdge;
  use crate::representation::payload::timetree::EdgeTimetree;
  use crate::representation::payload::timetree::NodeTimetree;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use bio::io::newick;
  use eyre::Report;
  use maplit::btreemap;
  use petgraph::visit::EdgeRef;
  use std::collections::BTreeMap;
  use std::io::Cursor;
  use treetime_graph::edge::{HasBranchLength, TimeLength};
  use treetime_graph::node::Named;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  #[test]
  fn test_create_branch_distributions_input_mode_sets_time_length() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let clock_rate = 0.001; // 0.001 subs/site/year

    create_branch_distributions_input_mode(&graph, clock_rate)?;

    // Verify each edge has time_length = branch_length / clock_rate
    for edge_ref in graph.get_edges() {
      let edge = edge_ref.read_arc().payload().read_arc();
      let branch_length = edge.branch_length();
      let time_length = edge.time_length();

      if let Some(bl) = branch_length {
        let expected_time = bl / clock_rate;
        let actual_time = time_length.expect("time_length should be set when branch_length exists");
        assert_ulps_eq!(actual_time, expected_time, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_input_mode_newick_output_uses_time_lengths() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let clock_rate = 0.001;

    create_branch_distributions_input_mode(&graph, clock_rate)?;

    // EdgeTimetree.nwk_weight() returns time_length, so Newick output should show time values
    let newick_output = nwk_write_str(&graph, &NwkWriteOptions::default())?;

    // Parse output with bio::io::newick (independent from our parser) to verify correctness
    let parsed = newick::read(Cursor::new(&newick_output)).expect("bio::newick should parse our output");

    // Build map of node name -> incoming edge weight (branch length to parent)
    // bio::newick stores edge weights as f32, node names as String
    let mut branch_lengths: BTreeMap<String, f64> = BTreeMap::new();
    for edge in parsed.g.edge_references() {
      let target_name = &parsed.g[edge.target()];
      if !target_name.is_empty() {
        branch_lengths.insert(target_name.clone(), *edge.weight() as f64);
      }
    }

    // Expected time lengths: branch_length / clock_rate
    // 0.003/0.001=3, 0.006/0.001=6, 0.009/0.001=9, 0.012/0.001=12
    let expected: BTreeMap<&str, f64> = btreemap! {
      "A" => 3.0,
      "B" => 6.0,
      "AB" => 9.0,
      "C" => 12.0,
    };

    for (name, expected_length) in expected {
      let actual_length = branch_lengths
        .get(name)
        .unwrap_or_else(|| panic!("Node '{name}' not found in parsed tree"));
      assert_ulps_eq!(*actual_length, expected_length, max_ulps = 4);
    }

    Ok(())
  }

  /// Gamma=2.0 means the branch evolves twice as fast, so the same number of
  /// substitutions corresponds to half the time duration.
  /// time = branch_length / (clock_rate * gamma)
  #[test]
  fn test_input_mode_gamma_scales_time_length() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.006)I:0.003)root;")?;
    let clock_rate = 0.001;

    // Set gamma=2.0 on edge to A
    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref.read_arc();
      let target = edge_read.target();
      let target_name = graph
        .get_node(target)
        .and_then(|n| n.read_arc().payload().read_arc().name().map(|s| s.as_ref().to_owned()));
      if target_name.as_deref() == Some("A") {
        edge_read.payload().write_arc().set_gamma(2.0);
      }
    }

    create_branch_distributions_input_mode(&graph, clock_rate)?;

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref.read_arc();
      let target = edge_read.target();
      let target_name = graph
        .get_node(target)
        .and_then(|n| n.read_arc().payload().read_arc().name().map(|s| s.as_ref().to_owned()));
      let payload = edge_read.payload().read_arc();

      match target_name.as_deref() {
        Some("A") => {
          // branch_length=0.006, gamma=2.0: time = 0.006 / (0.001 * 2.0) = 3.0
          // Newick parsing introduces tiny float error in branch length, use abs_diff
          let expected = 3.0;
          let actual = payload.time_length().expect("time_length should be set");
          assert_abs_diff_eq!(actual, expected, epsilon = 1e-7);
        },
        Some("I") => {
          // branch_length=0.003, gamma=1.0 (default): time = 0.003 / 0.001 = 3.0
          let expected = 3.0;
          let actual = payload.time_length().expect("time_length should be set");
          assert_abs_diff_eq!(actual, expected, epsilon = 1e-7);
        },
        _ => {},
      }
    }

    Ok(())
  }

  /// With gamma=1.0 (default), input mode produces identical results to the
  /// pre-gamma behavior: time = branch_length / clock_rate.
  #[test]
  fn test_input_mode_gamma_default_matches_no_gamma() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let clock_rate = 0.001;

    // All edges have default gamma=1.0
    create_branch_distributions_input_mode(&graph, clock_rate)?;

    for edge_ref in graph.get_edges() {
      let edge = edge_ref.read_arc().payload().read_arc();
      if let Some(bl) = edge.branch_length() {
        // With gamma=1.0, time = bl / clock_rate (same as without gamma)
        let expected = bl / clock_rate;
        let actual = edge.time_length().expect("time_length should be set");
        assert_ulps_eq!(actual, expected, max_ulps = 4);
      }
    }

    Ok(())
  }
}
