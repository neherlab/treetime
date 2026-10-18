#[cfg(test)]
mod tests {
  use crate::coalescent::mean_field::LineageExpectation;
  use crate::test_utils::{find_node_key_by_name, point_coalescent_prior};
  use crate::timetree::inference::result::NodePosterior;
  use eyre::Report;
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::pretty_assert_abs_diff_eq;

  // Two tips sampled at the same time t_s merge at the root. Under the Kingman coalescent the root time has density
  // (1/Tc) exp(-(t_s - t)/Tc) for t < t_s: two lineages, one pair, merger rate 1/Tc. The other node's expectation
  // (the tips) does not depend on the root, so the mean-field factor is exact here for any point time of the root.
  #[rustfmt::skip]
  #[rstest]
  #[case::root_far_back(   2.0, 1990.0)]
  #[case::root_near_tips(  2.0, 2009.0)]
  #[case::long_time_scale(25.0, 2000.0)]
  #[trace]
  fn test_mean_field_cherry_root_density_is_kingman(#[case] tc: f64, #[case] root_time: f64) -> Result<(), Report> {
    let tip = 2010.0;
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let prior = point_coalescent_prior(&graph, &names, &[("A", tip), ("B", tip), ("root", root_time)], tc)?;
    let root = find_node_key_by_name(&graph, &names, "root").expect("root must exist");

    let reference = 2005.0;
    for t in [1960.0, 1995.0, 2004.0, 2009.9] {
      let actual = prior.node_contribution(root, 2, t)? - prior.node_contribution(root, 2, reference)?;
      pretty_assert_abs_diff_eq!((reference - t) / tc, actual, epsilon = 1e-9);
    }
    Ok(())
  }

  // ((A,B)I,C)root with point times: the expectation is the lineage count k(t). With Tc = 1 the slope of an internal
  // node's negative log factor is 1 - k_{-v}(t), where k_{-v} counts the lineages without the node's own merger. The
  // steps of kbar are resolved to one cell of the expectation grid, 1e-3 years here.
  #[test]
  fn test_mean_field_internal_slope_counts_the_lineages_without_the_own_merger() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)I:1,C:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let times = [
      ("A", 2015.0),
      ("B", 2014.0),
      ("C", 2013.0),
      ("I", 2011.0),
      ("root", 2005.0),
    ];
    let prior = point_coalescent_prior(&graph, &names, &times, 1.0)?;
    let internal = find_node_key_by_name(&graph, &names, "I").expect("I must exist");
    let slope = |a: f64, b: f64| -> Result<f64, Report> {
      Ok((prior.node_contribution(internal, 2, b)? - prior.node_contribution(internal, 2, a)?) / (b - a))
    };

    // between the root and C, without I's merger: A, B and C
    pretty_assert_abs_diff_eq!(1.0 - 3.0, slope(2006.0, 2012.0)?, epsilon = 1e-3);
    // between C and B: A and B
    pretty_assert_abs_diff_eq!(1.0 - 2.0, slope(2013.2, 2013.8)?, epsilon = 1e-3);
    // between B and A: A only
    pretty_assert_abs_diff_eq!(1.0 - 1.0, slope(2014.2, 2014.8)?, epsilon = 1e-3);
    // before the root: two lineages remain once I's merger is removed
    pretty_assert_abs_diff_eq!(1.0 - 2.0, slope(1990.0, 2000.0)?, epsilon = 1e-3);
    Ok(())
  }

  // Point times give the statistics of the lineage counts: ((A,B)I,C)root with the root at 2005, I at 2011 and tips at
  // 2013, 2014, 2015 has k = 2, 3, 2, 1 on the four intervals, so int k (k - 1) / 2 dt = 6 + 6 + 1 + 0 and two mergers.
  // The steps of kbar are resolved to one cell of the expectation grid.
  #[test]
  fn test_mean_field_segment_statistics_of_point_times_are_those_of_the_lineage_counts() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:1,B:1)I:1,C:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let times = [
      ("A", 2015.0),
      ("B", 2014.0),
      ("C", 2013.0),
      ("I", 2011.0),
      ("root", 2005.0),
    ];
    let expectation = helpers::expectation(&graph, &names, &times, &[])?;

    let (pair_integral, mergers) = expectation.segment_statistics(&[2005.0, 2015.0]);

    pretty_assert_abs_diff_eq!(13.0, pair_integral[0], epsilon = 2e-2);
    pretty_assert_abs_diff_eq!(2.0, mergers[0], epsilon = 1e-12);
    Ok(())
  }

  // A cherry with tips at 2010 and the root uniform on [2000, 2004]: k = 2 between the root and the tips, so
  // E[int k (k - 1) / 2 dt] = E[2010 - t_root] = 8 with one merger, and Tc = I / M = 8. The variance of k enters: at a
  // time t inside the root's range E[k (k - 1)] / 2 = P(t_root < t), not kbar (kbar - 1) / 2. The tips' step is resolved
  // to half a cell of the expectation grid, 1.25e-3 here.
  #[test]
  fn test_mean_field_segment_statistics_integrate_over_the_root_posterior() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let flat_root = Distribution::function(Array1::linspace(2000.0, 2004.0, 401), Array1::zeros(401))?;
    let expectation = helpers::expectation(&graph, &names, &[("A", 2010.0), ("B", 2010.0)], &[("root", flat_root)])?;

    let (pair_integral, mergers) = expectation.segment_statistics(&[2000.0, 2010.0]);

    pretty_assert_abs_diff_eq!(8.0, pair_integral[0], epsilon = 5e-3);
    pretty_assert_abs_diff_eq!(1.0, mergers[0], epsilon = 1e-9);
    Ok(())
  }

  mod helpers {
    use super::*;

    /// Expectation from point times of the `points` nodes and the posteriors of the `distributions` nodes
    pub(super) fn expectation(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      points: &[(&str, f64)],
      distributions: &[(&str, Distribution<NegLog>)],
    ) -> Result<Arc<LineageExpectation>, Report> {
      let key = |name: &str| find_node_key_by_name(graph, names, name).expect("node must exist");
      let posterior = points
        .iter()
        .map(|&(name, time)| {
          let entry = NodePosterior {
            time: Some(time),
            ..NodePosterior::default()
          };
          (key(name), entry)
        })
        .chain(distributions.iter().map(|(name, distribution)| {
          let entry = NodePosterior {
            distribution: Some(Arc::new(distribution.clone())),
            ..NodePosterior::default()
          };
          (key(name), entry)
        }))
        .collect();
      let bad_branches = graph.get_nodes().map(|node| (node.key(), false)).collect();
      LineageExpectation::from_posterior(graph, &posterior, &bad_branches)
    }
  }
}
