#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::{empty_time_inference, find_node_key_by_name};
  use crate::timetree::inference::time_inference::{NodePosterior, TimeInference, likely_times};
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  #[test]
  fn test_time_inference_likely_times_names_node_with_nan_distribution() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let inference = helpers::inference_with_posterior(&graph, key, helpers::nan_distribution()?);

    let expected = format!(
      "When finding the most likely time of node {key}: \
       Cannot find the most likely time of a distribution function: its values contain NaN"
    );
    assert_error!(
      likely_times(&graph, &DateConstraints::default(), Some(&inference)),
      expected
    );
    assert_error!(inference.coalescent_node_times(), expected);
    Ok(())
  }

  #[test]
  fn test_time_inference_likely_times_prefer_the_date_constraint() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let inference = helpers::inference_with_posterior(&graph, key, Distribution::point(2010.0, 0.0));
    let mut constraints = DateConstraints::default();
    constraints
      .date_constraints
      .insert(key, Some(Arc::new(Distribution::point(2005.0, 0.0))));

    let with_constraint = likely_times(&graph, &constraints, Some(&inference))?;
    let without_constraint = likely_times(&graph, &DateConstraints::default(), Some(&inference))?;
    let before_inference = likely_times(&graph, &DateConstraints::default(), None)?;

    assert_eq!(Some(2005.0), with_constraint[&key]);
    assert_eq!(Some(2010.0), without_constraint[&key]);
    assert_eq!(None, before_inference[&key]);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn nan_distribution() -> Result<Distribution<NegLog>, Report> {
      Distribution::function(array![0.0, 1.0, 2.0], array![1.0, f64::NAN, 2.0])
    }

    pub(super) fn inference_with_posterior(
      graph: &Graph,
      key: GraphNodeKey,
      distribution: Distribution<NegLog>,
    ) -> TimeInference {
      let mut inference = empty_time_inference(graph);
      inference.posterior.insert(
        key,
        NodePosterior {
          distribution: Some(Arc::new(distribution)),
          ..NodePosterior::default()
        },
      );
      inference
    }
  }
}
