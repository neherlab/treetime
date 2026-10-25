#[cfg(test)]
mod tests {
  use crate::optimize::dispatch::optimize_edge_lengths;
  use crate::optimize::indel::EdgeIndels;
  use crate::optimize::params::BranchOptMethod;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::key_map::EdgeMap;

  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::missing( None,                "missing branch length")]
  #[case::negative(Some(-0.1),          "branch length must be finite and non-negative, got -0.1")]
  #[case::nan(     Some(f64::NAN),      "branch length must be finite and non-negative, got NaN")]
  #[case::infinite(Some(f64::INFINITY), "branch length must be finite and non-negative, got inf")]
  #[trace]
  fn test_dispatch_input_lengths_rejects_unusable_length_and_leaves_input_unchanged(
    #[case] branch_length: Option<f64>,
    #[case] reason: &str,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;".as_slice())?;
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let edge_key = graph.get_edges().last().expect("tree has edges").key();
    branch_lengths.insert(edge_key, branch_length);
    let input = format!("{branch_lengths:?}");

    let result = optimize_edge_lengths(
      &graph,
      10,
      &EdgeMap::new(),
      &EdgeIndels::Excluded,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    );

    let expected = match branch_length {
      None => format!("Cannot optimize edge {edge_key} with a {reason}"),
      Some(_) => format!("Cannot optimize edge {edge_key}: {reason}"),
    };
    assert_error!(result, expected);
    assert_eq!(input, format!("{branch_lengths:?}"));
    Ok(())
  }
}
