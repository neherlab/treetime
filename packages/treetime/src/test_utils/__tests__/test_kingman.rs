#[cfg(test)]
mod tests {
  use crate::test_utils::{
    NodeAge, find_node_key_by_name, kingman_neg_log_density, node_ages_from_dates, pair_exposure,
  };
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::key_map::NodeMap;

  use treetime_io::nwk::nwk_read;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::binary(
    "((A:1,B:1)N:1,C:1)root;",
    &[("A", 2010.0), ("B", 2009.0), ("C", 2012.0), ("N", 2005.0), ("root", 2000.0)],
    (1.0 + 4.0 * 3.0 + 5.0) / 3.0 + 2.0 * 3.0_f64.ln(),
  )]
  #[case::polytomy(
    "(A:1,B:1,C:1)root;",
    &[("A", 2010.0), ("B", 2009.0), ("C", 2012.0), ("root", 2000.0)],
    (1.0 + 9.0 * 3.0) / 3.0 + 2.0 * 3.0_f64.ln(),
  )]
  #[trace]
  fn test_kingman_neg_log_density_counts_a_polytomy_as_simultaneous_mergers(
    #[case] tree: &str,
    #[case] times: &[(&str, f64)],
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let parsed = nwk_read(tree.as_bytes())?;
    let dates: NodeMap<_> = times
      .iter()
      .map(|&(name, date)| {
        let key = find_node_key_by_name(&parsed.graph, &parsed.names, name).expect("node must exist");
        (key, date)
      })
      .collect();

    let actual = kingman_neg_log_density(&node_ages_from_dates(&parsed.graph, &dates), 3.0)?;

    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_kingman_pair_exposure_rejects_ages_that_leave_several_lineages() {
    let leaves = [
      NodeAge {
        age: 0.0,
        n_children: 0,
      },
      NodeAge {
        age: 1.0,
        n_children: 0,
      },
    ];

    let result = pair_exposure(&leaves);

    assert_error!(
      result,
      "Node ages leave 2 lineages after the oldest node; a tree ends with one"
    );
  }
}
