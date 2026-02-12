#[cfg(test)]
mod tests {
  use crate::commands::timetree::inference::runner::create_branch_distributions_input_mode;
  use crate::graph::edge::{HasBranchLength, TimeLength};
  use crate::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use crate::representation::edge_timetree::EdgeTimetree;
  use crate::representation::node_timetree::NodeTimetree;
  use approx::assert_ulps_eq;
  use bio::io::newick;
  use eyre::Report;
  use maplit::btreemap;
  use petgraph::visit::EdgeRef;
  use std::collections::BTreeMap;
  use std::io::Cursor;

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
        assert_ulps_eq!(time_length.unwrap_or(f64::NAN), expected_time, max_ulps = 4);
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
}
