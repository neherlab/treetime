#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::progress::NoopProgress;
  use crate::test_utils::{FullPosterior, full_posterior, parent_edge_key, point_date_constraints};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::forward_pass::{ForwardPassInputs, propagate_distributions_forward};
  use crate::timetree::inference::result::{BranchLikelihood, NodeBackward, log_marginal_likelihood};
  use crate::timetree::inference::runner::{EPS, GRID_POINTS};
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, DistributionFunction, NegLog, rewindow_to_mass};
  use treetime_graph::graph::Graph;
  use treetime_graph::key_map::{EdgeMap, NodeMap};
  use treetime_graph::node_names::NodeNames;
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, GridFn, MaxGridPoints, Side, SoftTailLaw};
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::ToF64Exact;

  const SITES: f64 = 754.0;
  const CLOCK_RATE: f64 = 0.0014;
  const TIP_DATE: f64 = 2023.0;
  const BRANCH_GRID_END: f64 = 10.0;

  #[rustfmt::skip]
  #[rstest]
  #[case::tips_20(  20, [300,  599, 1197])]
  #[case::tips_50(  50, [300,  599, 1197])]
  #[case::tips_200(200, [1197, 2393, 4785])]
  #[trace]
  fn test_zero_length_ladder_log_marginal_likelihood_converges_at_second_order_to_the_exact_value(
    #[case] n_tips: usize,
    #[case] branch_points: [usize; 3],
  ) -> Result<(), Report> {
    let expected = helpers::exact_log_z_from_independent_gaps(n_tips);

    let errors: Vec<f64> = branch_points
      .iter()
      .map(|&points| -> Result<f64, Report> {
        let ladder = helpers::Ladder::new(n_tips, points)?;
        let backward = ladder.backward(MaxGridPoints::default())?;
        let actual = log_marginal_likelihood(&backward, ladder.graph.root_key()?, None)?.expect("the root has a subtree");
        Ok((actual - expected).abs())
      })
      .try_collect()?;

    let reductions = [errors[0] / errors[1], errors[1] / errors[2]];
    assert!(
      reductions.iter().all(|&reduction| reduction >= 3.0),
      "errors {errors:?} fall by {reductions:?} per halving of the branch grid spacing"
    );
    Ok(())
  }

  #[test]
  fn test_zero_length_ladder_of_400_tips_dates_every_node_within_the_default_grid_point_limit() -> Result<(), Report> {
    let ladder = helpers::Ladder::new(400, GRID_POINTS)?;

    let backward = ladder.backward(MaxGridPoints::default())?;
    let posterior = ladder.forward(backward, MaxGridPoints::default())?;

    let undated = posterior.values().filter(|node| node.time.is_none()).count();
    assert_eq!(0, undated);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct Ladder {
      pub(super) graph: Graph,
      names: NodeNames,
      branches: EdgeMap<BranchLikelihood>,
      constraints: DateConstraints,
    }

    impl Ladder {
      pub(super) fn new(n_tips: usize, branch_points: usize) -> Result<Self, Report> {
        let newick = caterpillar_newick(n_tips);
        let parsed = nwk_read(newick.as_bytes())?;
        let names = parsed.names.clone();
        let graph = parsed.graph;
        let branch = zero_mutation_branch(branch_points)?;
        let branches = graph
          .get_nodes()
          .filter(|node| !node.is_root())
          .map(|node| {
            let branch = BranchLikelihood {
              distribution: Some(Arc::clone(&branch)),
              time_length: Some(0.0),
            };
            (parent_edge_key(&graph, node.key()), branch)
          })
          .collect();
        let tip_names = (0..n_tips).map(|i| format!("T{i}")).collect_vec();
        let dates = tip_names.iter().map(|name| (name.as_str(), TIP_DATE)).collect_vec();
        let constraints = point_date_constraints(&graph, &names, &dates);
        Ok(Self {
          graph,
          names,
          branches,
          constraints,
        })
      }

      pub(super) fn backward(&self, max_points: MaxGridPoints) -> Result<NodeMap<NodeBackward>, Report> {
        propagate_distributions_backward(
          &self.graph,
          &self.constraints,
          None,
          &self.bad_branches()?,
          &self.branches,
          max_points,
        )
      }

      pub(super) fn forward(
        &self,
        backward: NodeMap<NodeBackward>,
        max_points: MaxGridPoints,
      ) -> Result<NodeMap<FullPosterior>, Report> {
        let bad_branches = self.bad_branches()?;
        let inputs = ForwardPassInputs {
          constraints: &self.constraints,
          coalescent: None,
          bad_branches: &bad_branches,
          names: &self.names,
          branches: &self.branches,
          max_points,
        };
        propagate_distributions_forward(&self.graph, &inputs, backward, full_posterior, &NoopProgress)
      }

      fn bad_branches(&self) -> Result<NodeMap<bool>, Report> {
        derive_bad_branches(
          &self.graph,
          &self.constraints,
          &bad_leaves(&self.graph, &self.constraints, &BTreeSet::new()),
        )
      }
    }

    pub(super) fn exact_log_z_from_independent_gaps(n_tips: usize) -> f64 {
      let lambda = SITES * CLOCK_RATE;
      (0..n_tips - 1).map(|j| -(lambda * (j + 2).to_f64_exact()).ln()).sum()
    }

    fn caterpillar_newick(n_tips: usize) -> String {
      let inner = (1..n_tips).fold(String::from("T0:0"), |subtree, i| {
        format!("({subtree},T{i}:0)N{}:0", n_tips - 1 - i)
      });
      let root = inner
        .strip_suffix(":0")
        .expect("the caterpillar ends with a branch length");
      format!("{root};")
    }

    fn zero_mutation_branch(points: usize) -> Result<Arc<Distribution<NegLog>>, Report> {
      let lambda = SITES * CLOCK_RATE;
      let t = Array1::linspace(0.0, BRANCH_GRID_END, points);
      let f = GridFn::from_range_values((0.0, BRANCH_GRID_END), t.mapv(|t| lambda * t))?;
      let right = SoftTailLaw::fit(&f, Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
      let distribution = DistributionFunction::from_grid_fn(f)
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Linear(right))?;
      let distribution = rewindow_to_mass(
        &Distribution::Function(distribution),
        EPS,
        points,
        MaxGridPoints::default(),
      )?;
      Ok(Arc::new(distribution))
    }
  }
}
