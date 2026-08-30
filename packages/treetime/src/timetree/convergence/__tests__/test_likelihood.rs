#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::storage::dense::{DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
  use crate::partition::timetree::{GraphTimetree, PartitionTimetree, PartitionTimetreeRef};
  use crate::payload::timetree::NodeTimetree;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::convergence::likelihood::{
    compute_coalescent_log_lh, compute_positional_log_lh, compute_sequence_log_lh,
  };
  use crate::timetree::convergence::node_times::NodeTimeChange;
  use crate::timetree::convergence::optimizer::TimetreeOptimizer;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::LogLh;
  use treetime_utils::{o, pretty_assert_ulps_eq};

  #[test]
  fn test_likelihood_sequence_log_lh_sums_root_components() -> Result<(), Report> {
    let (graph, root_key) = helpers::single_root_graph()?;
    let partitions = [
      helpers::partition_with_root_log_lh(root_key, -2.0)?,
      helpers::partition_with_root_log_lh(root_key, -3.5)?,
    ];
    // Oracle: graph_log_lh() contract in packages/treetime/src/partition/traits.rs.
    let expected = -5.5;

    let actual = compute_sequence_log_lh(&graph, &partitions)
      .expect("sequence log-likelihood must be available")
      .value();

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_likelihood_sequence_log_lh_absent_without_partitions() -> Result<(), Report> {
    let (graph, _) = helpers::single_root_graph()?;

    let actual = compute_sequence_log_lh(&graph, &[]);

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_likelihood_positional_log_lh_sums_log_probabilities() -> Result<(), Report> {
    let graph = helpers::positional_graph()?;
    // Oracle: one edge with probability 0.25 contributes ln(0.25).
    let expected = 0.25_f64.ln();

    let actual = compute_positional_log_lh(&graph)
      .expect("positional log-likelihood must be available")
      .value();

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_likelihood_positional_log_lh_absent_without_distributions() -> Result<(), Report> {
    let graph: GraphTimetree = nwk_read_str("(child:0.1)root;")?;

    let actual = compute_positional_log_lh(&graph);

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_likelihood_coalescent_log_lh_matches_total_lh() -> Result<(), Report> {
    let graph = helpers::coalescent_graph()?;
    let tc = Distribution::constant(1.0);
    // Oracle: compute_coalescent_total_lh() is the coalescent model's whole-tree log-likelihood.
    let expected = compute_coalescent_total_lh(&graph, &tc)?.value();

    let actual = compute_coalescent_log_lh(&graph, Some(&tc))
      .expect("coalescent log-likelihood must be available")
      .value();

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_likelihood_coalescent_log_lh_absent_without_model() -> Result<(), Report> {
    let graph = helpers::coalescent_graph()?;

    let actual = compute_coalescent_log_lh(&graph, None);

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_likelihood_optimizer_total_sums_available_log_lh_components() -> Result<(), Report> {
    let graph = helpers::positional_graph()?;
    let root_key = find_node_key_by_name(&graph, "root").expect("root must exist");
    let partitions = [helpers::partition_with_root_log_lh(root_key, -2.0)?];
    let mut optimizer = TimetreeOptimizer::new(1, false);
    let expected = -2.0 + 0.25_f64.ln();

    assert!(optimizer.next_iter().is_some());
    optimizer.record(1, 0, NodeTimeChange::default(), &graph, &partitions, None)?;
    let actual = optimizer
      .trace()
      .first()
      .expect("one convergence metric must be recorded")
      .log_lh_total
      .expect("total log-likelihood must be available")
      .value();

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 10);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub fn single_root_graph() -> Result<(GraphTimetree, GraphNodeKey), Report> {
      let mut graph = GraphTimetree::new();
      let root_key = graph.add_node(NodeTimetree::default());
      graph.build()?;
      Ok((graph, root_key))
    }

    pub fn partition_with_root_log_lh(root_key: GraphNodeKey, log_lh: f64) -> Result<PartitionTimetreeRef, Report> {
      let mut partition = PartitionMarginalDense::new(0, jc69(JC69Params::default())?, Alphabet::default(), 1);
      partition.data.nodes.insert(
        root_key,
        DenseNodePartition {
          seq: DenseSeqInfo::default(),
          profile: DenseSeqDistribution::new(array![[1.0, 0.0, 0.0, 0.0]], LogLh::new(log_lh)),
        },
      );
      Ok(Arc::new(RwLock::new(PartitionTimetree::Dense(partition))))
    }

    pub fn positional_graph() -> Result<GraphTimetree, Report> {
      let graph: GraphTimetree = nwk_read_str("(child:0.1)root;")?;
      let root_key = find_node_key_by_name(&graph, "root").expect("root must exist");
      let child_key = find_node_key_by_name(&graph, "child").expect("child must exist");
      graph
        .get_node(root_key)
        .expect("root must exist")
        .read_arc()
        .payload()
        .write_arc()
        .time = Some(2000.0);
      graph
        .get_node(child_key)
        .expect("child must exist")
        .read_arc()
        .payload()
        .write_arc()
        .time = Some(2005.0);
      graph
        .get_edges()
        .into_iter()
        .next()
        .expect("one edge must exist")
        .read_arc()
        .payload()
        .write_arc()
        .branch_length_distribution = Some(Arc::new(Distribution::range((0.0, 10.0), 0.25)));
      Ok(graph)
    }

    pub fn coalescent_graph() -> Result<GraphTimetree, Report> {
      let dates = btreemap! {
        o!("root") => Some(DateConstraint::exact(2000.0)),
        o!("internal1") => Some(DateConstraint::exact(2005.0)),
        o!("leaf1") => Some(DateConstraint::exact(2010.0)),
        o!("leaf2") => Some(DateConstraint::exact(2010.0)),
        o!("leaf3") => Some(DateConstraint::exact(2012.0)),
      };
      let graph: GraphTimetree = nwk_read_str("((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;")?;
      load_date_constraints(&dates, &graph)?;
      Ok(graph)
    }
  }
}
