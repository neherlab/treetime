#[cfg(test)]
mod tests {
  use crate::o;
  use crate::reroot::params::{RerootMethod, RerootSpec};
  use eyre::Report;
  use helpers::{Target, reroot_repeatedly};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::pretty_assert_map_ulps_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::outgroup_of_a_two_child_root(  "((M:0.4,N:0.1)X:0.2,(C:0.3,D:0.5)Y:0.1)R;", Target::Outgroup(&["M"]), btreemap! { o!("M") => 0.2,      o!("X") => 0.2 })]
  #[case::middle_of_the_unrooted_branch( "(M:1,(A:1,B:1)S:3)R;",                       Target::Outgroup(&["M"]), btreemap! { o!("M") => 2.0,      o!("S") => 2.0 })]
  #[case::middle_inside_the_edge_band(   "(M:1.000002,(A:1,B:1)S:1)R;",                Target::Outgroup(&["M"]), btreemap! { o!("M") => 1.000001, o!("S") => 1.000001 })]
  #[case::below_a_three_child_root(      "(M:1,(A:1,B:1)S:2,C:3)R;",                   Target::Outgroup(&["M"]), btreemap! { o!("M") => 0.5,      o!("R") => 0.5 })]
  #[case::oldest_tip(                    "(M:1,(A:1,B:1)S:3)R;",                       Target::Spec(RerootSpec::Method(RerootMethod::Oldest)), btreemap! { o!("M") => 2.0, o!("S") => 2.0 })]
  #[trace]
  fn test_named_node_roots_at_the_middle_of_the_unrooted_branch_and_keeps_the_root_on_repeat(
    #[case] newick: &str,
    #[case] target: Target<'static>,
    #[case] expected: BTreeMap<String, f64>,
  ) -> Result<(), Report> {
    let once = reroot_repeatedly(newick, target.clone(), 1)?;
    let twice = reroot_repeatedly(newick, target, 2)?;

    pretty_assert_map_ulps_eq!(expected, once.root_children, max_ulps = 4);
    assert_eq!(once.newick, twice.newick, "a second reroot on the same target keeps the root");
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

    #[derive(Clone, Debug)]
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
