#[cfg(test)]
pub(crate) mod tests {
  use crate::o;
  use eyre::Report;
  use helpers::{midpoint, midpoint_error, midpoint_with_dated};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::{pretty_assert_map_abs_diff_eq, pretty_assert_map_ulps_eq};

  #[test]
  fn test_midpoint_star_roots_on_the_longest_branch() -> Result<(), Report> {
    let actual = midpoint("(A:1,B:2,C:6)R;", 1)?;

    let root_at_one_third_of_the_c_branch = btreemap! { o!("C") => 4.0, o!("R") => 2.0 };
    pretty_assert_map_ulps_eq!(root_at_one_third_of_the_c_branch, actual.root_children, max_ulps = 4);
    let expected = btreemap! { o!("A") => 3.0, o!("B") => 4.0, o!("C") => 4.0 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_midpoint_cherry_roots_at_equal_distances() -> Result<(), Report> {
    let actual = midpoint("(A:1,B:3)R;", 1)?;

    let expected = btreemap! { o!("A") => 2.0, o!("B") => 2.0 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_midpoint_inside_the_edge_end_band_is_placed_exactly() -> Result<(), Report> {
    let actual = midpoint("(A:1,B:2,C:2.00002)R;", 1)?;

    let root_1e_5_from_r = btreemap! { o!("C") => 2.00001, o!("R") => 0.00001 };
    pretty_assert_map_abs_diff_eq!(root_1e_5_from_r, actual.root_children, epsilon = 1e-15);
    let expected = btreemap! { o!("A") => 1.00001, o!("B") => 2.00001, o!("C") => 2.00001 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_midpoint_on_an_inner_node_roots_at_that_node() -> Result<(), Report> {
    let actual = midpoint("(A:1,(B:1,(C:1,D:1)Y:1)X:1)R;", 1)?;

    assert_eq!("X", actual.root_name);
    let expected = btreemap! { o!("A") => 2.0, o!("B") => 1.0, o!("C") => 2.0, o!("D") => 2.0 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::star(         "(A:1,B:2,C:6)R;")]
  #[case::inner_node(   "(A:1,(B:1,(C:1,D:1)Y:1)X:1)R;")]
  #[case::undated_stem( "((A:1,B:3)V:100)R;")]
  #[case::arms_equal_within_rounding("(T1:0.2888943450814206,(T0:0.001,T2:0.5308770876289646):0.20611932474743935)R;")]
  #[trace]
  fn test_midpoint_second_reroot_keeps_the_root(#[case] newick: &str) -> Result<(), Report> {
    let once = midpoint(newick, 1)?;

    let twice = midpoint(newick, 2)?;

    assert_eq!(
      (once.root_name, once.root_children, once.root_to_tip),
      (twice.root_name, twice.root_children, twice.root_to_tip)
    );
    Ok(())
  }

  #[test]
  fn test_midpoint_removes_an_undated_stem_root() -> Result<(), Report> {
    let actual = midpoint("((A:1,B:3)V:100)R;", 1)?;

    let expected = btreemap! { o!("A") => 2.0, o!("B") => 2.0 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_midpoint_keeps_a_dated_stem_root_as_a_tip() -> Result<(), Report> {
    let actual = midpoint_with_dated("((A:1,B:3)V:100)R;", &["R"], 1)?;

    let root_48_5_above_v = btreemap! { o!("R") => 51.5, o!("V") => 48.5 };
    pretty_assert_map_ulps_eq!(root_48_5_above_v, actual.root_children, max_ulps = 4);
    let expected = btreemap! { o!("A") => 49.5, o!("B") => 51.5, o!("R") => 51.5 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_midpoint_tie_of_the_longest_path_takes_the_first_in_postorder() -> Result<(), Report> {
    let actual = midpoint("(A:2,(B:2,C:2)X:0)R;", 1)?;

    assert_eq!("X", actual.root_name, "X closes a longest path before R in postorder");
    let expected = btreemap! { o!("A") => 2.0, o!("B") => 2.0, o!("C") => 2.0 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_midpoint_works_with_zero_length_branches() -> Result<(), Report> {
    let actual = midpoint("((A:0,B:2)X:0,C:2)R;", 1)?;

    let expected = btreemap! { o!("A") => 0.0, o!("B") => 2.0, o!("C") => 2.0 };
    pretty_assert_map_ulps_eq!(expected, actual.root_to_tip, max_ulps = 4);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::missing(  "(A:1,B)R;",    "R -> B: missing")]
  #[case::negative( "(A:1,B:-1)R;", "R -> B: -1")]
  #[trace]
  fn test_midpoint_rejects_invalid_branch_lengths(#[case] newick: &str, #[case] invalid_edge: &str) {
    let actual = midpoint_error(newick);

    assert_eq!(
      format!(
        "--reroot=midpoint places the root from branch lengths and requires every edge to have a finite, non-negative \
         branch length. Invalid edges:\n  {invalid_edge}\nGive the tree branch lengths, for example with `treetime \
         optimize`"
      ),
      actual
    );
  }

  pub(crate) mod helpers {
    use crate::branch_lengths::branch_length_or_zero;
    use crate::clock::divergence::root_to_node_divergences;
    use crate::progress::NoopProgress;
    use crate::reroot::date_free::{LengthsAdvice, reroot_date_free};
    use crate::reroot::params::{RerootMethod, RerootSpec};
    use crate::test_utils::find_node_key_by_name;
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::common_ancestor::common_ancestor;
    use treetime_graph::key_map::EdgeMap;

    use treetime_graph::graph::Graph;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};
    use treetime_utils::error::report_to_string;

    pub(crate) struct Rerooted {
      pub(crate) root_name: String,
      pub(crate) root_children: BTreeMap<String, f64>,
      pub(crate) root_to_tip: BTreeMap<String, f64>,
      pub(crate) diameter: f64,
      pub(crate) newick: String,
    }

    pub(crate) fn midpoint(newick: &str, count: usize) -> Result<Rerooted, Report> {
      midpoint_with_dated(newick, &[], count)
    }

    pub(super) fn midpoint_with_dated(newick: &str, dated: &[&str], count: usize) -> Result<Rerooted, Report> {
      let parsed = nwk_read(newick.as_bytes())?;
      let mut names = parsed.names;
      let mut graph = parsed.graph;
      let mut branch_lengths = parsed.branch_lengths;
      let dated: BTreeSet<_> = dated
        .iter()
        .map(|name| find_node_key_by_name(&graph, &names, name).expect("the dated node exists"))
        .collect();
      for _ in 0..count {
        reroot_date_free(
          &mut graph,
          &mut branch_lengths,
          &RerootSpec::Method(RerootMethod::Midpoint),
          &|key| dated.contains(&key),
          LengthsAdvice::OptimizeTheTree,
          &mut names,
          &NoopProgress,
        )?;
      }
      let schedule = GraphSchedule::new(&graph)?;
      let root_key = graph.get_exactly_one_root()?.key();
      let root_children = schedule
        .children(root_key)
        .iter()
        .map(|&(child_key, edge_key)| {
          (
            names[child_key].to_owned(),
            branch_length_or_zero(&branch_lengths, edge_key),
          )
        })
        .collect();
      let divergences = root_to_node_divergences(&graph, |edge_key| branch_length_or_zero(&branch_lengths, edge_key))?;
      let root_to_tip = schedule
        .leaves()
        .iter()
        .map(|&key| (names[key].to_owned(), divergences[&key]))
        .collect();
      let newick = nwk_write_str(
        &schedule.into_tree()?,
        &names,
        &branch_lengths,
        &NwkWriteOptions::default(),
        &NwkNodeComments::new(),
      )?;
      Ok(Rerooted {
        root_name: names[root_key].to_owned(),
        root_children,
        root_to_tip,
        diameter: diameter(&graph, &branch_lengths)?,
        newick,
      })
    }

    pub(crate) fn diameter(graph: &Graph, branch_lengths: &EdgeMap<Option<f64>>) -> Result<f64, Report> {
      let divergences = root_to_node_divergences(graph, |edge_key| branch_length_or_zero(branch_lengths, edge_key))?;
      let leaves = graph.get_leaves().map(|leaf| leaf.key()).collect_vec();
      leaves
        .iter()
        .tuple_combinations()
        .map(|(&a, &b)| Ok(divergences[&a] + divergences[&b] - 2.0 * divergences[&common_ancestor(graph, &[a, b])?]))
        .try_fold(0.0_f64, |longest, distance: Result<f64, Report>| {
          Ok(longest.max(distance?))
        })
    }

    pub(super) fn midpoint_error(newick: &str) -> String {
      let Err(report) = midpoint(newick, 1) else {
        panic!("midpoint rooting must reject the branch lengths of {newick}");
      };
      report_to_string(&report)
    }
  }
}
