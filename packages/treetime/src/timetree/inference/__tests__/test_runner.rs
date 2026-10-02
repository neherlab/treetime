#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::pretty_assert_ulps_eq;
  use crate::timetree::inference::result::BranchLikelihood;
  use crate::timetree::inference::runner::create_branch_distributions_input_mode;
  use crate::timetree::optimization::relaxed_clock::unit_gammas;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

  #[test]
  fn test_create_branch_distributions_input_mode_sets_time_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.003,B:0.006)AB:0.009,C:0.012)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let clock_rate = 0.001;

    let branches = create_branch_distributions_input_mode(&graph, &branch_lengths, &unit_gammas(&graph), clock_rate);

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref;
      let key = edge_read.key();
      let branch_length = branch_lengths.get(&key).copied().flatten();
      let time_length = branches[&key].time_length;

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

    let branches = create_branch_distributions_input_mode(&graph, &branch_lengths, &unit_gammas(&graph), clock_rate);

    let time_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        (key, branches[&key].time_length)
      })
      .collect();
    let actual = nwk_write_str(&graph, &names, &time_lengths, &NwkWriteOptions::default())?;

    let expected = "((A:3,B:6)AB:9,C:12)root;";
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_input_mode_gamma_scales_time_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.006)I:0.003)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let clock_rate = 0.001;

    let gammas = graph
      .get_edges()
      .map(|edge| {
        let target_name = names[&edge.target()].as_deref();
        let gamma = if target_name == Some("A") { 2.0 } else { 1.0 };
        (edge.key(), gamma)
      })
      .collect();

    let branches = create_branch_distributions_input_mode(&graph, &branch_lengths, &gammas, clock_rate);

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref;
      let key = edge_read.key();
      let target = edge_read.target();
      let target_name = graph
        .get_node(target)
        .and_then(|n| names.get(&n.key()).cloned().flatten());
      let time_length = branches[&key].time_length;

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

    let branches = create_branch_distributions_input_mode(&graph, &branch_lengths, &unit_gammas(&graph), clock_rate);

    for edge_ref in graph.get_edges() {
      let edge_read = edge_ref;
      let key = edge_read.key();
      if let Some(bl) = branch_lengths.get(&key).copied().flatten() {
        let expected = bl / clock_rate;
        let actual = branches[&key].time_length.expect("time_length should be set");
        pretty_assert_ulps_eq!(actual, expected, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_input_mode_edge_without_branch_length_has_no_branch_likelihood() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A)root;")?;
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let edge_key = graph
      .get_edges()
      .collect::<Vec<_>>()
      .pop()
      .expect("tree must contain one edge")
      .key();
    branch_lengths.insert(edge_key, None);

    let branches = create_branch_distributions_input_mode(&graph, &branch_lengths, &unit_gammas(&graph), 0.001);

    let expected = BranchLikelihood {
      distribution: None,
      time_length: None,
    };
    assert_eq!(expected, branches[&edge_key]);
    Ok(())
  }
}
