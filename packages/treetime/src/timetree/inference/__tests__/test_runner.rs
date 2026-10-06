#[cfg(test)]
mod tests {
  use crate::o;
  use crate::progress::NoopProgress;
  use crate::test_utils::{RecordingLog, find_node_key_by_name, parent_edge_key};
  use crate::timetree::inference::result::{BranchLikelihood, NodeTimes};
  use crate::timetree::inference::runner::{
    CLOCK_BRANCH_LENGTH_DAMPING, blended_clock_branch_lengths, create_branch_distributions_input_mode,
    explain_grid_point_limit,
  };
  use crate::timetree::optimization::relaxed_clock::unit_gammas;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::tree_view::TreeView;
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};
  use treetime_utils::error::report_to_string;
  use treetime_utils::make_report;

  const DYADIC_CLOCK_RATE: f64 = 0.5;

  const BLEND_TREE: &str = "((A:1,B:1)AB:1,C:1,D:1)root;";

  #[test]
  fn test_runner_grid_point_limit_error_names_the_fixed_clock_rate() {
    let exceeded = MaxGridPoints::default()
      .point_count(15_583_861.0, (2019.0, 2022.4), 2.2e-7)
      .unwrap_err()
      .wrap_err("When sending the time message backward along edge 3");
    let expected = "Time inference needs a grid of 15583861 points, more than the limit of 1000000 \
      (--max-grid-points). The fixed --clock-rate 0.5 may be far from the rate the data supports, which makes \
      branch-time distributions very narrow. Check --clock-rate, or raise --max-grid-points when enough memory is \
      available.: When sending the time message backward along edge 3: A grid over [2019, 2022.4] with spacing \
      2.2e-7 needs 15583861 points, more than the limit of 1000000";
    assert_eq!(
      expected,
      report_to_string(&explain_grid_point_limit(exceeded, 0.5, true))
    );
  }

  #[test]
  fn test_runner_grid_point_limit_error_names_the_estimated_clock_rate() {
    let exceeded = MaxGridPoints::default()
      .point_count(f64::INFINITY, (0.0, 1.0), 1e-320)
      .unwrap_err();
    let expected = "Time inference needs a grid with more points than the limit of 1000000 (--max-grid-points). \
      The estimated clock rate 0.00123 may be far from the true rate, which makes branch-time distributions very \
      narrow; errors in the input dates and weak temporal signal cause such estimates. Check the input dates, set \
      --clock-rate, or raise --max-grid-points when enough memory is available.: A grid over [0, 1] with spacing \
      1.0e-320 needs more points than the limit of 1000000";
    assert_eq!(
      expected,
      report_to_string(&explain_grid_point_limit(exceeded, 1.234e-3, false))
    );
  }

  #[test]
  fn test_runner_grid_point_limit_explanation_leaves_other_errors_unchanged() {
    let other = make_report!("Clock rate is negative");
    assert_eq!(
      "Clock rate is negative",
      report_to_string(&explain_grid_point_limit(other, 0.5, true))
    );
  }

  #[test]
  fn test_create_branch_distributions_input_mode_gives_each_edge_a_point_at_its_time_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:1.5,B:3.0)AB:4.5,C:6.0)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let branches = create_branch_distributions_input_mode(
      &graph,
      &nwk_parsed.branch_lengths,
      &unit_gammas(&graph),
      DYADIC_CLOCK_RATE,
    );

    let expected = btreemap! {
      o!("A") => helpers::point_branch(3.0),
      o!("AB") => helpers::point_branch(9.0),
      o!("B") => helpers::point_branch(6.0),
      o!("C") => helpers::point_branch(12.0),
    };
    assert_eq!(expected, helpers::by_target_name(&graph, &names, &branches));
    Ok(())
  }

  #[test]
  fn test_input_mode_newick_output_uses_time_lengths() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.003,B:0.006)AB:0.009,C:0.012)root;".as_slice())?;
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
    let actual = nwk_write_str(
      &TreeView::new(&graph)?,
      &names,
      &time_lengths,
      &NwkWriteOptions::default(),
      &NwkNodeComments::new(),
    )?;

    let expected = "((A:3,B:6)AB:9,C:12)root;";
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_input_mode_gamma_scales_time_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0)I:1.5)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let gammas = btreemap! {
      parent_edge_key(&graph, key("A")) => 2.0,
      parent_edge_key(&graph, key("I")) => 1.0,
    };

    let branches =
      create_branch_distributions_input_mode(&graph, &nwk_parsed.branch_lengths, &gammas, DYADIC_CLOCK_RATE);

    let expected = btreemap! {
      o!("A") => helpers::point_branch(3.0),
      o!("I") => helpers::point_branch(3.0),
    };
    assert_eq!(expected, helpers::by_target_name(&graph, &names, &branches));
    Ok(())
  }

  #[test]
  fn test_input_mode_edge_without_branch_length_has_no_branch_likelihood() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A)root;".as_slice())?;
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

  #[test]
  fn test_blended_clock_branch_lengths_blends_scales_clamps_and_drops_edges() -> Result<(), Report> {
    let fixture = helpers::BlendFixture::new()?;

    let actual = fixture.blend(&NoopProgress);

    let expected = btreemap! {
      o!("A") => 6.0,
      o!("AB") => 1.5,
      o!("B") => 0.0,
      o!("D") => 7.0,
    };
    assert_eq!(expected, fixture.by_target_name(&actual));
    Ok(())
  }

  #[test]
  fn test_blended_clock_branch_lengths_warns_about_inverted_branches() -> Result<(), Report> {
    let fixture = helpers::BlendFixture::new()?;
    let log = RecordingLog::default();

    fixture.blend(&log);

    let expected = vec![
      "Timetree: 1 branch(es) run backwards in time, i.e. the child is dated before its parent. Their clock \
       branch lengths were committed as zero. This is expected only where an observed leaf date or an exact \
       internal-node date conflicts with the fitted clock, since the forward pass clamps the other internal nodes \
       to their parent but leaves exact dates as given."
        .to_owned(),
    ];
    assert_eq!(expected, log.warnings());
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::progress::LogSink;

    pub(super) fn point_branch(time_length: f64) -> BranchLikelihood {
      BranchLikelihood {
        distribution: Some(Arc::new(Distribution::point(time_length, 0.0))),
        time_length: Some(time_length),
      }
    }

    pub(super) fn by_target_name<V: Clone>(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      values: &BTreeMap<GraphEdgeKey, V>,
    ) -> BTreeMap<String, V> {
      graph
        .get_edges()
        .filter_map(|edge| {
          let name = names[&edge.target()].clone().expect("every fixture node is named");
          let value = values.get(&edge.key())?.clone();
          Some((name, value))
        })
        .collect()
    }

    pub(super) struct BlendFixture {
      graph: Graph,
      names: BTreeMap<GraphNodeKey, Option<String>>,
      previous: BTreeMap<GraphEdgeKey, f64>,
      node_times: NodeTimes,
      gammas: BTreeMap<GraphEdgeKey, f64>,
    }

    impl BlendFixture {
      pub(super) fn new() -> Result<Self, Report> {
        let nwk_parsed = nwk_read(BLEND_TREE.as_bytes())?;
        let names = nwk_parsed.names();
        let graph = nwk_parsed.graph;
        let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
        let edge = |name: &str| parent_edge_key(&graph, key(name));
        let removed_edge = GraphEdgeKey(graph.get_edges().count() + 10);
        let previous = btreemap! {
          edge("AB") => 1.0,
          edge("D") => 7.0,
          removed_edge => 3.0,
        };
        let node_times = btreemap! {
          key("root") => Some(2000.0),
          key("AB") => Some(2004.0),
          key("A") => Some(2010.0),
          key("B") => Some(2002.0),
          key("C") => None,
          key("D") => None,
        };
        let mut gammas = unit_gammas(&graph);
        gammas.insert(edge("A"), 2.0);
        Ok(Self {
          graph,
          names,
          previous,
          node_times,
          gammas,
        })
      }

      pub(super) fn blend(&self, log: &dyn LogSink) -> BTreeMap<GraphEdgeKey, f64> {
        blended_clock_branch_lengths(
          &self.graph,
          DYADIC_CLOCK_RATE,
          CLOCK_BRANCH_LENGTH_DAMPING,
          &self.previous,
          &self.node_times,
          &self.gammas,
          log,
        )
      }

      pub(super) fn by_target_name(&self, values: &BTreeMap<GraphEdgeKey, f64>) -> BTreeMap<String, f64> {
        assert!(
          values.keys().all(|key| self.graph.get_edge(*key).is_some()),
          "every committed length must belong to a current graph edge"
        );
        by_target_name(&self.graph, &self.names, values)
      }
    }
  }
}
