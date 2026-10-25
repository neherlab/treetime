#[cfg(test)]
mod tests {
  use crate::o;
  use crate::reroot::params::{RerootMethod, RerootSpec};
  use eyre::Report;
  use helpers::{Target, reroot_repeatedly};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_utils::pretty_assert_map_ulps_eq;

  #[test]
  fn test_named_node_second_reroot_on_an_outgroup_keeps_the_root() -> Result<(), Report> {
    let newick = "((M:0.4,N:0.1)X:0.2,(C:0.3,D:0.5)Y:0.1)R;";

    let once = reroot_repeatedly(newick, Target::Outgroup(&["M"]), 1)?;
    let twice = reroot_repeatedly(newick, Target::Outgroup(&["M"]), 2)?;

    let expected = btreemap! { o!("M") => 0.2, o!("X") => 0.2 };
    pretty_assert_map_ulps_eq!(expected, once.root_children, max_ulps = 4);
    assert_eq!(
      once.newick, twice.newick,
      "a second reroot on the same outgroup keeps the root"
    );
    Ok(())
  }

  #[test]
  fn test_named_node_child_of_a_two_child_root_roots_at_the_middle_of_its_unrooted_branch() -> Result<(), Report> {
    let newick = "(M:1,(A:1,B:1)S:3)R;";

    let actual = reroot_repeatedly(newick, Target::Outgroup(&["M"]), 1)?;

    let expected = btreemap! { o!("M") => 2.0, o!("S") => 2.0 };
    pretty_assert_map_ulps_eq!(expected, actual.root_children, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_named_node_middle_inside_the_edge_end_band_is_placed_exactly() -> Result<(), Report> {
    let newick = "(M:1.000002,(A:1,B:1)S:1)R;";

    let actual = reroot_repeatedly(newick, Target::Outgroup(&["M"]), 1)?;

    let expected = btreemap! { o!("M") => 1.000001, o!("S") => 1.000001 };
    pretty_assert_map_ulps_eq!(expected, actual.root_children, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_named_node_second_reroot_below_a_three_child_root_keeps_the_root() -> Result<(), Report> {
    let newick = "(M:1,(A:1,B:1)S:2,C:3)R;";

    let once = reroot_repeatedly(newick, Target::Outgroup(&["M"]), 1)?;
    let twice = reroot_repeatedly(newick, Target::Outgroup(&["M"]), 2)?;

    let expected = btreemap! { o!("M") => 0.5, o!("R") => 0.5 };
    pretty_assert_map_ulps_eq!(expected, once.root_children, max_ulps = 4);
    assert_eq!(
      once.newick, twice.newick,
      "a second reroot on the same outgroup keeps the root"
    );
    Ok(())
  }

  #[test]
  fn test_named_node_oldest_tip_roots_at_the_middle_of_its_unrooted_branch() -> Result<(), Report> {
    let newick = "(M:1,(A:1,B:1)S:3)R;";

    let once = reroot_repeatedly(newick, Target::Spec(RerootSpec::Method(RerootMethod::Oldest)), 1)?;
    let twice = reroot_repeatedly(newick, Target::Spec(RerootSpec::Method(RerootMethod::Oldest)), 2)?;

    let expected = btreemap! { o!("M") => 2.0, o!("S") => 2.0 };
    pretty_assert_map_ulps_eq!(expected, once.root_children, max_ulps = 4);
    assert_eq!(
      once.newick, twice.newick,
      "a second reroot on the oldest tip keeps the root"
    );
    Ok(())
  }

  mod helpers {
    use crate::clock::clock_state::ClockInputs;
    use crate::progress::NoopProgress;
    use crate::reroot::driver::reroot_at_best_root;
    use crate::reroot::params::{RerootParams, RerootSpec};
    use crate::test_utils::find_node_key_by_name;
    use eyre::Report;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::key_map::EdgeMap;
    use treetime_graph::key_map::NodeMap;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};

    pub(super) enum Target<'a> {
      Outgroup(&'a [&'a str]),
      Spec(RerootSpec),
    }

    pub(super) struct RerootedTree {
      pub(super) root_children: BTreeMap<String, f64>,
      pub(super) newick: String,
    }

    pub(super) fn reroot_repeatedly(newick: &str, target: Target<'_>, count: usize) -> Result<RerootedTree, Report> {
      let parsed = nwk_read(newick.as_bytes())?;
      let mut names = parsed.names;
      let mut graph = parsed.graph;
      let mut branch_lengths = parsed.branch_lengths;
      let spec = match target {
        Target::Outgroup(tips) => RerootSpec::Tips(
          tips
            .iter()
            .map(|tip| find_node_key_by_name(&graph, &names, tip).expect("the tip exists"))
            .collect(),
        ),
        Target::Spec(spec) => spec,
      };
      let dates: NodeMap<_> = [
        ("M", 2000.0),
        ("A", 2010.0),
        ("B", 2011.0),
        ("C", 2012.0),
        ("N", 2013.0),
      ]
      .into_iter()
      .filter_map(|(name, date)| Some((find_node_key_by_name(&graph, &names, name)?, Some(date))))
      .collect();
      let mut inputs = ClockInputs::from_times(&graph, &dates, &EdgeMap::new());
      for _ in 0..count {
        reroot_at_best_root(
          &mut graph,
          &mut branch_lengths,
          &mut inputs,
          &BTreeSet::new(),
          &RerootParams::new(spec.clone(), false),
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
          let length = branch_lengths[&edge_key].expect("every edge has a length");
          (names[child_key].to_owned(), length)
        })
        .collect();
      let newick = nwk_write_str(
        &schedule.into_tree()?,
        &names,
        &branch_lengths,
        &NwkWriteOptions::default(),
        &NwkNodeComments::new(),
      )?;
      Ok(RerootedTree { root_children, newick })
    }
  }
}
