#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::coalescent::node_time::CoalescentNodeTimes;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::shared::update::MarginalEdges;
  use crate::partition::storage::dense::{DenseNodeState, DenseSeqDistribution, DenseSeqInfo};
  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::NodeSeqInput;
  use crate::test_utils::find_node_key_by_name;
  use crate::test_utils::{constraint_coalescent_node_times, empty_time_inference};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::convergence::likelihood::{
    compute_coalescent_log_lh, compute_positional_log_lh, compute_sequence_log_lh,
  };
  use crate::timetree::convergence::node_times::NodeTimeChange;
  use crate::timetree::convergence::optimizer::TimetreeOptimizer;
  use crate::timetree::inference::time_inference::{BranchLikelihood, TimeInference};
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{LogLh, seq};
  use treetime_utils::{o, pretty_assert_ulps_eq};

  #[test]
  fn test_likelihood_sequence_log_lh_is_the_root_log_lh_of_the_partition() -> Result<(), Report> {
    let (graph, root_key) = helpers::single_root_graph()?;
    let branch_model = BranchModel::Marginal(helpers::partition_with_root_log_lh(&graph, root_key, -3.5)?);
    let expected = -3.5;

    let actual = compute_sequence_log_lh(&graph, &branch_model)
      .expect("sequence log-likelihood must be available")
      .value();

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_likelihood_sequence_log_lh_absent_in_input_mode() -> Result<(), Report> {
    let (graph, _) = helpers::single_root_graph()?;

    let actual = compute_sequence_log_lh(&graph, &BranchModel::Input);

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_likelihood_positional_log_lh_sums_log_probabilities() -> Result<(), Report> {
    let (graph, names) = helpers::positional_graph()?;
    let expected = 0.25_f64.ln();

    let state = helpers::positional_state(&graph, &names);
    let actual = compute_positional_log_lh(&graph, &state)
      .expect("positional log-likelihood must be available")
      .value();

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_likelihood_positional_log_lh_absent_without_distributions() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(child:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;

    let state = empty_time_inference(&graph);
    let actual = compute_positional_log_lh(&graph, &state);

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_likelihood_coalescent_log_lh_matches_total_lh() -> Result<(), Report> {
    let (graph, constraints) = helpers::coalescent_graph()?;
    let tc = Distribution::constant(1.0);
    let node_times = helpers::coalescent_node_times(&graph, &constraints);
    let expected = compute_coalescent_total_lh(&graph, &tc, &node_times, &BTreeMap::new(), &NoopProgress)?.value();

    let actual = compute_coalescent_log_lh(&graph, Some(&tc), &node_times, &BTreeMap::new(), &NoopProgress)
      .expect("coalescent log-likelihood must be available")
      .value();

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 10);
    Ok(())
  }

  #[test]
  fn test_likelihood_coalescent_log_lh_absent_without_model() -> Result<(), Report> {
    let (graph, constraints) = helpers::coalescent_graph()?;
    let node_times = helpers::coalescent_node_times(&graph, &constraints);

    let actual = compute_coalescent_log_lh(&graph, None, &node_times, &BTreeMap::new(), &NoopProgress);

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_likelihood_optimizer_total_sums_available_log_lh_components() -> Result<(), Report> {
    let (graph, names) = helpers::positional_graph()?;
    let root_key = find_node_key_by_name(&graph, &names, "root").expect("root must exist");
    let branch_model = BranchModel::Marginal(helpers::partition_with_root_log_lh(&graph, root_key, -2.0)?);
    let mut optimizer = TimetreeOptimizer::new(1, false);
    let expected = -2.0 + 0.25_f64.ln();
    let state = helpers::positional_state(&graph, &names);

    assert!(optimizer.next_iter(&NoopProgress).is_some());
    optimizer.record(
      1,
      0,
      NodeTimeChange::default(),
      &graph,
      &branch_model,
      &state,
      None,
      helpers::fixed_clock(),
      &BTreeMap::new(),
      &NoopProgress,
    )?;
    let actual = optimizer
      .trace
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
    use crate::timetree::convergence::metrics::IterationClock;

    pub(super) const fn fixed_clock() -> IterationClock {
      IterationClock {
        clock_rate: 1e-3,
        r_squared: None,
      }
    }

    pub(super) fn single_root_graph() -> Result<(Graph, GraphNodeKey), Report> {
      let mut graph = Graph::new();
      let root_key = graph.add_node();
      graph.build()?;
      Ok((graph, root_key))
    }

    pub(super) fn partition_with_root_log_lh(
      graph: &Graph,
      root_key: GraphNodeKey,
      log_lh: f64,
    ) -> Result<PartitionTimetree, Report> {
      let alphabet = Alphabet::default();
      let fill = alphabet.char(0);
      let node_inputs = graph
        .get_leaves()
        .map(|leaf| {
          (
            leaf.key(),
            NodeSeqInput {
              name: None,
              seq: Some(seq![fill]),
            },
          )
        })
        .collect();
      let partition = PartitionMarginalDense::new(0, alphabet, graph, &node_inputs)?;
      let node_states = btreemap! {
        root_key => DenseNodeState {
          seq: DenseSeqInfo::default(),
          profile: DenseSeqDistribution::new(array![[1.0, 0.0, 0.0, 0.0]], LogLh::new(log_lh)),
        },
      };
      Ok(PartitionTimetree::Dense(DenseReconstruction {
        partition,
        gtr: jc69(JC69Params::default())?,
        node_states,
        edges: MarginalEdges::default(),
      }))
    }

    pub(super) fn positional_graph() -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>), Report> {
      let nwk_parsed = nwk_read_str("(child:0.1)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let graph: Graph = graph;
      Ok((graph, names))
    }

    pub(super) fn positional_state(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> TimeInference {
      let root_key = find_node_key_by_name(graph, names, "root").expect("root must exist");
      let child_key = find_node_key_by_name(graph, names, "child").expect("child must exist");
      let mut state = empty_time_inference(graph);
      state
        .posterior
        .get_mut(&root_key)
        .expect("root must have a posterior")
        .time = Some(2000.0);
      state
        .posterior
        .get_mut(&child_key)
        .expect("child must have a posterior")
        .time = Some(2005.0);
      let edge_key = graph.get_edges().next().expect("one edge must exist").key();
      let branch = BranchLikelihood {
        distribution: Some(Arc::new(Distribution::range((0.0, 10.0), 0.25))),
        time_length: None,
      };
      state.branches.insert(edge_key, branch);
      state
    }

    pub(super) fn coalescent_graph() -> Result<(Graph, DateConstraints), Report> {
      let dates = btreemap! {
        o!("root") => Some(DateConstraint::exact(2000.0)),
        o!("internal1") => Some(DateConstraint::exact(2005.0)),
        o!("leaf1") => Some(DateConstraint::exact(2010.0)),
        o!("leaf2") => Some(DateConstraint::exact(2010.0)),
        o!("leaf3") => Some(DateConstraint::exact(2012.0)),
      };
      let nwk_parsed = nwk_read_str("((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let graph: Graph = graph;
      let constraints = load_date_constraints(&dates, &graph, &names, &NoopProgress)?;
      Ok((graph, constraints))
    }

    pub(super) fn coalescent_node_times(graph: &Graph, constraints: &DateConstraints) -> CoalescentNodeTimes {
      constraint_coalescent_node_times(graph, constraints).unwrap()
    }
  }
}
