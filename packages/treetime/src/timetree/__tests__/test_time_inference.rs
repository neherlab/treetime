#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::{empty_time_inference, find_node_key_by_name, point_date_constraints};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::result::{NodePosterior, TimeInference, given_times, likely_times};
  use crate::timetree::optimization::polytomy::resolve::{PolytomyResolution, resolve_polytomies};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;
  use treetime_utils::sync::random::get_random_number_generator;

  const POLYTOMY_MUTATION_RATE: f64 = 0.1;

  const POLYTOMY_MERGER_RATE: f64 = 0.15;

  const POLYTOMY_SEED: u64 = 11;

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
    assert_error!(likely_times(&graph, &DateConstraints::default(), &inference), expected);
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
      .by_node
      .insert(key, Some(Arc::new(Distribution::point(2005.0, 0.0))));

    let with_constraint = likely_times(&graph, &constraints, &inference)?;
    let without_constraint = likely_times(&graph, &DateConstraints::default(), &inference)?;
    let before_inference = given_times(&graph, &DateConstraints::default())?;

    assert_eq!(Some(2005.0), with_constraint[&key]);
    assert_eq!(Some(2010.0), without_constraint[&key]);
    assert_eq!(None, before_inference[&key]);
    Ok(())
  }

  #[test]
  fn test_time_inference_node_times_reads_the_committed_time_of_every_node() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let mut inference = empty_time_inference(&graph);
    for (name, time) in [
      ("A", Some(2010.0)),
      ("B", Some(2011.0)),
      ("I", Some(2005.0)),
      ("root", None),
    ] {
      inference.posterior.insert(
        key(name),
        NodePosterior {
          distribution: Some(Arc::new(Distribution::point(1900.0, 0.0))),
          time,
          contradicted: false,
        },
      );
    }

    let expected = btreemap! {
      key("A") => Some(2010.0),
      key("B") => Some(2011.0),
      key("I") => Some(2005.0),
      key("root") => None,
    };
    assert_eq!(expected, inference.node_times());
    Ok(())
  }

  #[test]
  fn test_time_inference_coalescent_node_times_pair_the_committed_time_with_the_distribution_peak() -> Result<(), Report>
  {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let mut inference = empty_time_inference(&graph);
    inference.posterior.insert(
      key("I"),
      NodePosterior {
        distribution: Some(Arc::new(Distribution::range((2004.0, 2006.0), 0.0))),
        time: Some(2006.5),
        contradicted: false,
      },
    );
    inference.posterior.insert(
      key("A"),
      NodePosterior {
        distribution: None,
        time: Some(2010.0),
        contradicted: false,
      },
    );
    inference.bad_branches.insert(key("B"), true);

    let actual = inference
      .coalescent_node_times()?
      .into_iter()
      .map(|(key, entry)| (key, (entry.time, entry.time_dist_likely, entry.bad_branch)))
      .collect::<BTreeMap<_, _>>();

    let expected = btreemap! {
      key("A") => (Some(2010.0), None, false),
      key("B") => (None, None, true),
      key("I") => (Some(2006.5), Some(2005.0), false),
      key("root") => (None, None, false),
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_time_inference_likely_times_cover_every_node() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let mut constraints = point_date_constraints(&graph, &names, &[("A", 2010.0)]);
    constraints
      .by_node
      .insert(key("B"), Some(Arc::new(Distribution::range((2008.0, 2012.0), 0.0))));
    let mut inference = empty_time_inference(&graph);
    for (name, peak) in [("A", 2000.0), ("I", 2005.0)] {
      inference.posterior.insert(
        key(name),
        NodePosterior {
          distribution: Some(Arc::new(Distribution::point(peak, 0.0))),
          ..NodePosterior::default()
        },
      );
    }

    let actual = likely_times(&graph, &constraints, &inference)?;

    let expected = btreemap! {
      key("A") => Some(2010.0),
      key("B") => Some(2010.0),
      key("I") => Some(2005.0),
      key("root") => None,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_time_inference_likely_times_take_the_posterior_of_a_polytomy_merger() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2,C:0.15)ABC:0.05)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let [a, b, c, abc, root] = ["A", "B", "C", "ABC", "root"].map(key);
    let times = btreemap! {
      a => Some(2020.0),
      b => Some(2015.0),
      c => Some(2018.0),
      abc => Some(1990.0),
      root => Some(1980.0),
    };
    let constraints = point_date_constraints(&graph, &names, &[("A", 2020.0), ("B", 2015.0), ("C", 2018.0)]);
    let branch_lengths = graph.get_edges().map(|edge| (edge.key(), Some(0.0))).collect();
    let PolytomyResolution {
      graph, merger_times, ..
    } = resolve_polytomies(
      graph,
      branch_lengths,
      &BranchModel::Input,
      POLYTOMY_MUTATION_RATE,
      0,
      &PiecewiseConstantFn::new(array![], array![POLYTOMY_MERGER_RATE]),
      &mut get_random_number_generator(POLYTOMY_SEED),
      &times,
    )?;
    let merger = *merger_times.keys().exactly_one().expect("the seed resolves one merger");
    let mut inference = empty_time_inference(&graph);
    inference.posterior.insert(
      merger,
      NodePosterior {
        distribution: Some(Arc::new(Distribution::point(2001.0, 0.0))),
        ..NodePosterior::default()
      },
    );

    let actual = likely_times(&graph, &constraints, &inference)?;

    let expected = btreemap! {
      a => Some(2020.0),
      b => Some(2015.0),
      c => Some(2018.0),
      abc => None,
      root => None,
      merger => Some(2001.0),
    };
    assert_eq!(expected, actual);
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
