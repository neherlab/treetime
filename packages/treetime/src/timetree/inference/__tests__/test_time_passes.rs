#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::progress::NoopProgress;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, point_date_constraints};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::timetree::inference::result::{BranchLikelihood, NodePosterior, TimeBackward};
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::collections::BTreeSet;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((B:1,C:1,(U1:1,U2:1)N:1)P:1,A:1)root;";

  const BRANCH_MEAN: f64 = 2.0;

  const BRANCH_SD: f64 = 0.5;

  const BRANCH_GRID_POINTS: usize = 201;

  #[test]
  fn test_time_passes_node_with_only_undated_leaves_takes_the_parent_message() -> Result<(), Report> {
    let mut fixture = helpers::Fixture::new()?;
    fixture.set_point_branch("N", BRANCH_MEAN);
    let constraints = fixture.dated(&[("A", 2011.0), ("B", 2010.0), ("C", 2010.5), ("P", 2008.0)]);

    let (bad_branches, backward, posterior) = fixture.run(&constraints)?;

    assert!(bad_branches[&fixture.key("N")], "N has only undated leaves below it");
    assert_eq!(None, backward.subtree[&fixture.key("N")]);
    let expected = NodePosterior {
      distribution: Some(Arc::new(Distribution::point(2008.0 + BRANCH_MEAN, 0.0))),
      likely_time: Some(2008.0 + BRANCH_MEAN),
      time: Some(2008.0 + BRANCH_MEAN),
      contradicted: false,
    };
    assert_eq!(expected, posterior[&fixture.key("N")]);
    Ok(())
  }

  #[test]
  fn test_time_passes_internal_node_with_its_own_date_keeps_the_date() -> Result<(), Report> {
    let fixture = helpers::Fixture::new()?;
    let constraints = fixture.dated(&[("A", 2011.0), ("B", 2010.0), ("C", 2010.5), ("N", 2008.0)]);

    let (bad_branches, backward, posterior) = fixture.run(&constraints)?;

    assert!(
      !bad_branches[&fixture.key("N")],
      "a node with its own date is never bad"
    );
    assert!(
      backward.messages[&fixture.edge("N")].is_some(),
      "the date of N must reach its parent as evidence"
    );
    let expected = NodePosterior {
      distribution: Some(Arc::new(Distribution::point(2008.0, 0.0).normalize()?)),
      likely_time: Some(2008.0),
      time: Some(2008.0),
      contradicted: false,
    };
    assert_eq!(expected, posterior[&fixture.key("N")]);
    Ok(())
  }

  #[test]
  fn test_time_passes_date_of_an_internal_node_pulls_its_parent_toward_it() -> Result<(), Report> {
    let fixture = helpers::Fixture::new()?;
    let leaf_dates = [("A", 2011.0), ("B", 2010.0), ("C", 2010.5)];
    let with_date_of_n = fixture.dated(&[leaf_dates.as_slice(), &[("N", 2008.0)]].concat());
    let without_date_of_n = fixture.dated(&leaf_dates);
    let mean_of_b_and_c_messages = f64::midpoint(2010.0, 2010.5) - BRANCH_MEAN;
    let mean_of_b_c_and_n_messages = (2010.0 + 2010.5 + 2008.0) / 3.0 - BRANCH_MEAN;

    let (_, with_backward, _) = fixture.run(&with_date_of_n)?;
    let (_, without_backward, _) = fixture.run(&without_date_of_n)?;

    let peak_of = |backward: &TimeBackward| -> Result<f64, Report> {
      let subtree = backward.subtree[&fixture.key("P")]
        .as_deref()
        .expect("P must have a subtree distribution");
      Ok(subtree.likely_time()?.expect("P must have a likely time"))
    };
    let with = peak_of(&with_backward)?;
    let without = peak_of(&without_backward)?;
    assert!(
      with < without,
      "the earlier date of N must pull P earlier: {with} vs {without}"
    );
    assert!(
      (with - mean_of_b_c_and_n_messages).abs() < (with - mean_of_b_and_c_messages).abs(),
      "with the date of N, P must peak near {mean_of_b_c_and_n_messages}, got {with}"
    );
    assert!(
      (without - mean_of_b_and_c_messages).abs() < (without - mean_of_b_c_and_n_messages).abs(),
      "without the date of N, P must peak near {mean_of_b_and_c_messages}, got {without}"
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct Fixture {
      graph: Graph,
      names: BTreeMap<GraphNodeKey, Option<String>>,
      branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
    }

    impl Fixture {
      pub(super) fn new() -> Result<Self, Report> {
        let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
        let names = nwk_parsed.names();
        let graph = nwk_parsed.graph;
        let branches = graph
          .get_edges()
          .map(|edge| {
            let branch = BranchLikelihood {
              distribution: Some(Arc::new(branch_distribution(BRANCH_MEAN)?)),
              time_length: Some(BRANCH_MEAN),
            };
            Ok((edge.key(), branch))
          })
          .collect::<Result<_, Report>>()?;
        Ok(Self { graph, names, branches })
      }

      pub(super) fn key(&self, name: &str) -> GraphNodeKey {
        find_node_key_by_name(&self.graph, &self.names, name).expect("fixture node must exist")
      }

      pub(super) fn edge(&self, name: &str) -> GraphEdgeKey {
        parent_edge_key(&self.graph, self.key(name))
      }

      pub(super) fn set_point_branch(&mut self, name: &str, length: f64) {
        let branch = BranchLikelihood {
          distribution: Some(Arc::new(Distribution::point(length, 0.0))),
          time_length: Some(length),
        };
        self.branches.insert(self.edge(name), branch);
      }

      pub(super) fn dated(&self, dates: &[(&str, f64)]) -> DateConstraints {
        point_date_constraints(&self.graph, &self.names, dates)
      }

      pub(super) fn run(
        &self,
        constraints: &DateConstraints,
      ) -> Result<
        (
          BTreeMap<GraphNodeKey, bool>,
          TimeBackward,
          BTreeMap<GraphNodeKey, NodePosterior>,
        ),
        Report,
      > {
        let bad_branches = derive_bad_branches(
          &self.graph,
          constraints,
          &bad_leaves(&self.graph, constraints, &BTreeSet::new()),
        )?;
        let backward = propagate_distributions_backward(&self.graph, constraints, None, &bad_branches, &self.branches)?;
        let posterior = propagate_distributions_forward(
          &self.graph,
          constraints,
          &self.names,
          &self.branches,
          &backward,
          &NoopProgress,
        )?;
        Ok((bad_branches, backward, posterior))
      }
    }

    fn branch_distribution(mean: f64) -> Result<Distribution<NegLog>, Report> {
      let t = Array1::linspace(0.0, 2.0 * mean, BRANCH_GRID_POINTS);
      let y = t.mapv(|t: f64| 0.5 * ((t - mean) / BRANCH_SD).powi(2));
      Distribution::function(t, y)
    }
  }
}
