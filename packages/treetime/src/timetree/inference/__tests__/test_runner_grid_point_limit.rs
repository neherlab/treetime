#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, unknown_branches};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::grid_request::{
    GridOperand, GridOperation, GridOrigin, GridRequest, OperandRole, WrapGridRequest,
  };
  use crate::timetree::inference::result::{BranchLikelihood, NodeBackward};
  use crate::timetree::inference::runner::explain_grid_point_limit;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::key_map::NodeMap;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::node_names::NodeNames;
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::error::report_to_string;
  use treetime_utils::{assert_error, make_report};

  const CLOCK_RATE: f64 = 0.5;

  const LEAF_TREE: &str = "(A:1)root;";

  const CHAIN_TREE: &str = "((A:1)B:1)root;";

  #[test]
  fn test_runner_grid_point_limit_names_the_date_constraint_that_sets_the_spacing() -> Result<(), Report> {
    let fixture = helpers::Fixture::new(LEAF_TREE)?;
    let a = fixture.key("A");
    let dates = vec![(a, helpers::gaussian(2000.0, 2000.25, 17)?)];
    let branches = vec![(a, helpers::gaussian(0.0, 16.0, 33)?)];

    let result = fixture.backward_pass(dates, branches, true);

    let expected = format!(
      "Time inference needs a grid of 1025 points, more than the limit of 1000 (--max-grid-points) when sending \
       the time message of node 'A' backward to its parent. The grid spans [-16, 0] at spacing 0.015625: the span is \
       that of the time distribution of the branch above node 'A', and the spacing is that of the time distribution \
       of node 'A' given its subtree, which takes it from the date constraint of node 'A'. The operands are the time \
       distribution of node 'A' given its subtree over [2000, 2000.25] at spacing 0.015625, and the time \
       distribution of the branch above node 'A' over [-16, 0] at spacing 0.5. Raise --max-grid-points when enough \
       memory is available.: When sending the time message of node {a} backward to its parent: A grid over [-16, 0] \
       with spacing 0.015625 needs 1025 points, more than the limit of 1000"
    );
    assert_error!(result, expected);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::fixed_rate(    true,  "The fixed --clock-rate 0.5 may be far above the rate the data supports, which makes branch-time distributions very narrow. Check --clock-rate")]
  #[case::estimated_rate(false, "The estimated clock rate 0.5 may be far above the true rate, which makes branch-time distributions very narrow; errors in the input dates and weak temporal signal cause such estimates. Check the input dates, set --clock-rate")]
  #[trace]
  fn test_runner_grid_point_limit_names_the_clock_rate_when_a_branch_sets_the_spacing(
    #[case] clock_rate_fixed: bool,
    #[case] advice: &str,
  ) -> Result<(), Report> {
    let fixture = helpers::Fixture::new(LEAF_TREE)?;
    let a = fixture.key("A");
    let dates = vec![(a, helpers::gaussian(2000.0, 2016.0, 33)?)];
    let branches = vec![(a, helpers::gaussian(0.0, 0.25, 17)?)];

    let result = fixture.backward_pass(dates, branches, clock_rate_fixed);

    let expected = format!(
      "Time inference needs a grid of 1025 points, more than the limit of 1000 (--max-grid-points) when sending \
       the time message of node 'A' backward to its parent. The grid spans [2000, 2016] at spacing 0.015625: the \
       span is that of the time distribution of node 'A' given its subtree, which takes it from the date constraint \
       of node 'A', and the spacing is that of the time distribution of the branch above node 'A'. The operands are \
       the time distribution of node 'A' given its subtree over [2000, 2016] at spacing 0.5, and the time \
       distribution of the branch above node 'A' over [-0.25, 0] at spacing 0.015625. {advice}, or raise \
       --max-grid-points when enough memory is available.: When sending the time message of node {a} backward to \
       its parent: A grid over [2000, 2016] with spacing 0.015625 needs 1025 points, more than the limit of 1000"
    );
    assert_error!(result, expected);
    Ok(())
  }

  #[test]
  fn test_runner_grid_point_limit_traces_a_node_spacing_to_the_branch_below_it_and_the_span_to_its_own_branch()
  -> Result<(), Report> {
    let fixture = helpers::Fixture::new(CHAIN_TREE)?;
    let a = fixture.key("A");
    let b = fixture.key("B");
    let dates = vec![(a, Distribution::point(2010.0, 0.0))];
    let branches = vec![
      (a, helpers::gaussian(0.0, 0.25, 17)?),
      (b, helpers::gaussian(0.0, 16.0, 33)?),
    ];

    let result = fixture.backward_pass(dates, branches, true);

    let expected = format!(
      "Time inference needs a grid of 1025 points, more than the limit of 1000 (--max-grid-points) when sending \
       the time message of node 'B' backward to its parent. The grid spans [-16, 0] at spacing 0.015625: the span is \
       that of the time distribution of the branch above node 'B', and the spacing is that of the time distribution \
       of node 'B' given its subtree, which takes it from the time distribution of the branch above node 'A'. The \
       operands are the time distribution of node 'B' given its subtree over [2009.75, 2010] at spacing 0.015625, \
       and the time distribution of the branch above node 'B' over [-16, 0] at spacing 0.5. The fixed --clock-rate \
       0.5 may be far from the rate the data supports, which makes branch-time distributions very narrow or very \
       wide. Check --clock-rate, or raise --max-grid-points when enough memory is available.: When sending the time message of \
       node {b} backward to its parent: A grid over [-16, 0] with spacing 0.015625 needs 1025 points, more than the \
       limit of 1000"
    );
    assert_error!(result, expected);
    Ok(())
  }

  #[test]
  fn test_runner_grid_point_limit_names_a_wide_branch_that_sets_the_span_of_an_overlap() -> Result<(), Report> {
    let fixture = helpers::Fixture::new(LEAF_TREE)?;
    let a = fixture.key("A");
    let root = fixture.key("root");
    let message = helpers::gaussian(2000.0, 2010.0, 11)?;
    let subtree = helpers::gaussian(2005.0, 2020.0, 16)?;

    let result = MaxGridPoints::new(1000)?
      .point_count(5001.0, (2005.0, 2010.0), 0.001)
      .wrap_grid_request(|| {
        let operands = [
          GridOperand::of(OperandRole::ForwardMessage(a), &message, GridOrigin::branch(a)),
          GridOperand::of(OperandRole::Subtree(a), &subtree, GridOrigin::date_constraint(root)),
        ];
        GridRequest::new(GridOperation::Posterior(a), operands)
      })
      .map_err(|report| explain_grid_point_limit(report, &fixture.names, CLOCK_RATE, false));

    let expected = format!(
      "Time inference needs a grid of 5001 points, more than the limit of 1000 (--max-grid-points) when multiplying \
       the forward time message of node 'A' with its time distribution given its subtree. The grid spans [2005, \
       2010] at spacing 0.001: the span is the overlap of the operands, and no operand has this spacing. The \
       operands are the forward time message to node 'A' over [2000, 2010] at spacing 1, and the time distribution \
       of node 'A' given its subtree over [2005, 2020] at spacing 1. The estimated clock rate 0.5 may be far below \
       the true rate, which makes branch-time distributions very wide; errors in the input dates and weak temporal \
       signal cause such estimates. Check the input dates, set --clock-rate, or raise --max-grid-points when enough \
       memory is available.: When multiplying the forward time message of node {a} with its time distribution given \
       its subtree: A grid over [2005, 2010] with spacing 0.001 needs 5001 points, more than the limit of 1000"
    );
    assert_error!(result, expected);
    Ok(())
  }

  #[test]
  fn test_runner_grid_point_limit_without_a_request_gives_no_clock_rate_advice() -> Result<(), Report> {
    let fixture = helpers::Fixture::new(LEAF_TREE)?;

    let result = MaxGridPoints::default()
      .point_count(f64::INFINITY, (0.0, 1.0), 1e-320)
      .map_err(|report| explain_grid_point_limit(report, &fixture.names, CLOCK_RATE, false));

    let expected = "Time inference needs a grid with more points than the limit of 1000000 (--max-grid-points). Raise \
      --max-grid-points when enough memory is available.: A grid over [0, 1] with spacing 1.0e-320 needs more \
      points than the limit of 1000000";
    assert_error!(result, expected);
    Ok(())
  }

  #[test]
  fn test_runner_grid_point_limit_explanation_leaves_other_errors_unchanged() -> Result<(), Report> {
    let fixture = helpers::Fixture::new(LEAF_TREE)?;
    let other = make_report!("Clock rate is negative");
    assert_eq!(
      "Clock rate is negative",
      report_to_string(&explain_grid_point_limit(other, &fixture.names, CLOCK_RATE, true))
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct Fixture {
      pub graph: Graph,
      pub names: NodeNames,
    }

    impl Fixture {
      pub(super) fn new(tree: &str) -> Result<Self, Report> {
        let nwk_parsed = nwk_read(tree.as_bytes())?;
        Ok(Self {
          names: nwk_parsed.names,
          graph: nwk_parsed.graph,
        })
      }

      pub(super) fn key(&self, name: &str) -> GraphNodeKey {
        find_node_key_by_name(&self.graph, &self.names, name).expect("fixture node must exist")
      }

      pub(super) fn backward_pass(
        &self,
        dates: Vec<(GraphNodeKey, Distribution<NegLog>)>,
        branches_above: Vec<(GraphNodeKey, Distribution<NegLog>)>,
        clock_rate_fixed: bool,
      ) -> Result<NodeMap<NodeBackward>, Report> {
        let mut constraints = DateConstraints::default();
        for (node, date) in dates {
          constraints.by_node.insert(node, Some(Arc::new(date)));
        }
        let mut branches = unknown_branches(&self.graph);
        for (node, branch) in branches_above {
          let branch = BranchLikelihood {
            distribution: Some(Arc::new(branch)),
            time_length: None,
          };
          branches.insert(parent_edge_key(&self.graph, node), branch);
        }
        let bad_branches = self.graph.get_nodes().map(|node| (node.key(), false)).collect();
        propagate_distributions_backward(
          &self.graph,
          &constraints,
          None,
          &bad_branches,
          &branches,
          MaxGridPoints::new(1000)?,
        )
        .map_err(|report| explain_grid_point_limit(report, &self.names, CLOCK_RATE, clock_rate_fixed))
      }
    }

    pub(super) fn gaussian(start: f64, end: f64, n_points: usize) -> Result<Distribution<NegLog>, Report> {
      let x = Array1::linspace(start, end, n_points);
      let mid = f64::midpoint(start, end);
      let y = x.mapv(|t| (t - mid).powi(2));
      Distribution::function(x, y)
    }
  }
}
