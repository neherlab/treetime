#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::progress::NoopProgress;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::timetree::inference::result::{BranchLikelihood, NodePosterior, TimeBackward};
  use crate::timetree::inference::runner::{EPS, GRID_POINTS};
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::collections::BTreeSet;
  use std::sync::Arc;
  use treetime_distribution::{
    Distribution, NegLog, convolve_across_edge, distribution_division, distribution_multiplication,
    distribution_product,
  };
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::Side;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((B:1,C:1,(U1:1,U2:1)N:1)P:1,A:1)root;";

  #[test]
  fn test_time_passes_node_with_only_undated_leaves_takes_the_parent_message() -> Result<(), Report> {
    let fixture = helpers::Fixture::new()?;
    let constraints = fixture.dated(&[("A", 2011.0), ("B", 2010.0), ("C", 2010.5)]);

    let (bad_branches, backward, posterior) = fixture.run(&constraints)?;

    assert!(bad_branches[&fixture.key("N")], "N has only undated leaves below it");
    assert_eq!(None, backward.subtree[&fixture.key("N")]);

    let expected = fixture.expected_posterior_of_n_from_its_parent()?;
    let actual = posterior[&fixture.key("N")]
      .distribution
      .as_deref()
      .expect("N must get a posterior from its parent");
    assert_eq!(&expected, actual);

    let parent_time = posterior[&fixture.key("P")].time.expect("P must be dated");
    let expected_time = expected
      .likely_time()?
      .map(|time| time.max(parent_time))
      .expect("the parent message must have a likely time");
    assert_eq!(Some(expected_time), posterior[&fixture.key("N")].time);
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
      time: Some(2008.0),
      contradicted: false,
    };
    assert_eq!(expected, posterior[&fixture.key("N")]);
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
              distribution: Some(Arc::new(branch_distribution(2.0)?)),
              time_length: Some(2.0),
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
        let key = self.key(name);
        self
          .graph
          .get_edges()
          .find(|edge| edge.target() == key)
          .expect("fixture node must have a parent edge")
          .key()
      }

      pub(super) fn dated(&self, dates: &[(&str, f64)]) -> DateConstraints {
        let date_constraints = dates
          .iter()
          .map(|(name, date)| (self.key(name), Some(Arc::new(Distribution::point(*date, 0.0)))))
          .collect();
        DateConstraints {
          by_node: date_constraints,
        }
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

      pub(super) fn expected_posterior_of_n_from_its_parent(&self) -> Result<Distribution<NegLog>, Report> {
        let branch = |name: &str| -> &Distribution<NegLog> {
          self.branches[&self.edge(name)]
            .distribution
            .as_deref()
            .expect("every fixture edge has a branch likelihood")
        };
        let message_up = |subtree: &Distribution<NegLog>, name: &str| -> Result<Distribution<NegLog>, Report> {
          convolve_across_edge(subtree, &branch(name).negate()?, Side::Left, EPS, GRID_POINTS)
        };
        let leaf = |date: f64| Distribution::point(date, 0.0).normalize();

        let msg_a = message_up(&leaf(2011.0)?, "A")?;
        let msg_b = message_up(&leaf(2010.0)?, "B")?;
        let msg_c = message_up(&leaf(2010.5)?, "C")?;
        let subtree_p = distribution_product(&[&msg_b, &msg_c])?.normalize()?;
        let msg_p = message_up(&subtree_p, "P")?;
        let posterior_root = distribution_product(&[&msg_p, &msg_a])?.normalize()?;

        let cavity_p = distribution_division(&posterior_root, &msg_p)?;
        let from_root = convolve_across_edge(&cavity_p, branch("P"), Side::Right, EPS, GRID_POINTS)?;
        let posterior_p = distribution_multiplication(&from_root, &subtree_p)?.normalize()?;

        convolve_across_edge(&posterior_p, branch("N"), Side::Right, EPS, GRID_POINTS)
      }
    }

    fn branch_distribution(mean: f64) -> Result<Distribution<NegLog>, Report> {
      let t = Array1::linspace(0.0, 2.0 * mean, 201);
      let y = t.mapv(|t: f64| 0.5 * ((t - mean) / 0.5).powi(2));
      Distribution::function(t, y)
    }
  }
}
