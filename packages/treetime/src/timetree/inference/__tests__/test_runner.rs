#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::pretty_assert_ulps_eq;
  use crate::timetree::inference::runner::create_branch_distributions_input_mode;
  use crate::timetree::timetree_state::TimetreeState;
  use approx::assert_abs_diff_eq;
  use bio::io::newick;
  use eyre::Report;
  use maplit::btreemap;
  use petgraph::visit::EdgeRef;
  use std::collections::BTreeMap;
  use std::io::Cursor;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  #[test]
  fn test_create_branch_distributions_input_mode_sets_time_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let clock_rate = 0.001;

    let mut state = TimetreeState::new(&graph);
    create_branch_distributions_input_mode(&graph, &branch_lengths, clock_rate, &mut state)?;

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref;
      let key = edge_read.key();
      let branch_length = branch_lengths.get(&key).copied().flatten();
      let time_length = state.edge(key).time_length;

      if let Some(bl) = branch_length {
        let expected_time = bl / clock_rate;
        let actual_time = time_length.expect("time_length should be set when branch_length exists");
        pretty_assert_ulps_eq!(actual_time, expected_time, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_input_mode_newick_output_uses_time_lengths() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let clock_rate = 0.001;

    let mut state = TimetreeState::new(&graph);
    create_branch_distributions_input_mode(&graph, &branch_lengths, clock_rate, &mut state)?;

    let time_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        (key, state.edge(key).time_length)
      })
      .collect();
    let newick_output = nwk_write_str(&graph, &names, &time_lengths, &NwkWriteOptions::default())?;

    let parsed = newick::read(Cursor::new(&newick_output)).expect("bio::newick should parse our output");

    let mut branch_lengths: BTreeMap<String, f64> = BTreeMap::new();
    for edge in parsed.g.edge_references() {
      let target_name = &parsed.g[edge.target()];
      if !target_name.is_empty() {
        branch_lengths.insert(target_name.clone(), *edge.weight() as f64);
      }
    }

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
      pretty_assert_ulps_eq!(*actual_length, expected_length, max_ulps = 4);
    }

    Ok(())
  }

  #[test]
  fn test_input_mode_gamma_scales_time_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.006)I:0.003)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let clock_rate = 0.001;

    let mut state = TimetreeState::new(&graph);

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref;
      let key = edge_read.key();
      let target = edge_read.target();
      let target_name = graph
        .get_node(target)
        .and_then(|n| names.get(&n.key()).cloned().flatten());
      if target_name.as_deref() == Some("A") {
        state.edge_mut(key).gamma = 2.0;
      }
    }

    create_branch_distributions_input_mode(&graph, &branch_lengths, clock_rate, &mut state)?;

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref;
      let key = edge_read.key();
      let target = edge_read.target();
      let target_name = graph
        .get_node(target)
        .and_then(|n| names.get(&n.key()).cloned().flatten());
      let time_length = state.edge(key).time_length;

      match target_name.as_deref() {
        Some("A") => {
          let expected = 3.0;
          let actual = time_length.expect("time_length should be set");
          assert_abs_diff_eq!(actual, expected, epsilon = 1e-7);
        },
        Some("I") => {
          let expected = 3.0;
          let actual = time_length.expect("time_length should be set");
          assert_abs_diff_eq!(actual, expected, epsilon = 1e-7);
        },
        _ => {},
      }
    }

    Ok(())
  }

  #[test]
  fn test_input_mode_gamma_default_matches_no_gamma() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let clock_rate = 0.001;

    let mut state = TimetreeState::new(&graph);
    create_branch_distributions_input_mode(&graph, &branch_lengths, clock_rate, &mut state)?;

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref;
      let key = edge_read.key();
      if let Some(bl) = branch_lengths.get(&key).copied().flatten() {
        let expected = bl / clock_rate;
        let actual = state.edge(key).time_length.expect("time_length should be set");
        pretty_assert_ulps_eq!(actual, expected, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_input_mode_uses_time_length_when_branch_length_is_absent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let edge = graph
      .get_edges()
      .collect::<Vec<_>>()
      .pop()
      .expect("tree must contain one edge");
    let edge_key = edge.key();
    branch_lengths.insert(edge_key, None);

    let mut state = TimetreeState::new(&graph);
    state.edge_mut(edge_key).time_length = Some(7.5);
    create_branch_distributions_input_mode(&graph, &branch_lengths, 0.001, &mut state)?;

    assert_eq!(Some(7.5), state.edge(edge_key).time_length);
    assert_eq!(
      Some(7.5),
      state
        .edge(edge_key)
        .branch_length_distribution
        .as_ref()
        .and_then(|distribution| distribution.likely_time().unwrap())
    );
    Ok(())
  }
}
