#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::time_inference::{BranchLikelihood, TimeBackward};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  #[test]
  fn test_backward_pass_computes_internal_node_time() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:2.5)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_key, 2013.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_key, 2.5);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let time_dist = node_time_distribution(&backward, internal_key)
      .expect("internal node should have time distribution after backward pass");
    let likely_time = time_dist
      .likely_time()?
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2010.5, max_ulps = 4);
    assert!(likely_time < 2013.0, "Parent should be older than child");

    Ok(())
  }

  #[test]
  fn test_backward_pass_nan_leaf_likelihood_error_reaches_caller() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:2.5)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let mut inputs = BackwardInputs::new(&graph);
    set_date_constraint(&mut inputs.constraints, leaf_key, Distribution::point(2013.0, f64::NAN));
    set_edge_branch_dist(&graph, &mut inputs, leaf_key, 2.5);

    let result = run_backward_pass(&graph, &inputs, None);

    assert_error!(
      result,
      format!(
        "When normalizing the time distribution of node {leaf_key}: Cannot normalize a distribution point: its peak negative log-likelihood is NaN"
      )
    );
    Ok(())
  }

  #[test]
  fn test_backward_pass_multiplies_child_messages() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_a_key, 2015.0);
    set_leaf_time(&mut inputs, leaf_b_key, 2014.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_b_key, 2.0);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let time_dist =
      node_time_distribution(&backward, internal_key).expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()?
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_backward_pass_preserves_leaf_time_distribution_with_coalescent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");

    let date_a = 2015.0;
    let date_b = 2014.0;

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_a_key, date_a);
    set_leaf_time(&mut inputs, leaf_b_key, date_b);
    set_edge_branch_dist(&graph, &mut inputs, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_b_key, 2.0);

    let coalescent_model = coalescent_model(0.01)?;

    let first = run_backward_pass(&graph, &inputs, Some(&coalescent_model))?;
    let backward = run_backward_pass(&graph, &inputs, Some(&coalescent_model))?;
    assert_eq!(first, backward);

    {
      let time_dist = node_time_distribution(&backward, leaf_a_key).expect("leaf A should have time distribution");
      let expected = Distribution::point(date_a, 0.0);
      assert_eq!(&expected, time_dist.as_ref());
    }

    {
      let time_dist = node_time_distribution(&backward, leaf_b_key).expect("leaf B should have time distribution");
      let expected = Distribution::point(date_b, 0.0);
      assert_eq!(&expected, time_dist.as_ref());
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_preserves_internal_time_with_strong_coalescent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal I not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_a_key, 2015.0);
    set_leaf_time(&mut inputs, leaf_b_key, 2014.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_b_key, 2.0);

    let coalescent_model = coalescent_model(1e-6)?;

    let backward = run_backward_pass(&graph, &inputs, Some(&coalescent_model))?;

    let actual =
      node_time_distribution(&backward, internal_key).and_then(|distribution| distribution.likely_time().unwrap());
    let expected = Some(2012.0);
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_backward_pass_leaf_subtree_distribution_is_its_date_constraint() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let constraint = Distribution::range((2014.0, 2015.0), 0.0);
    let mut inputs = BackwardInputs::new(&graph);
    set_date_constraint(&mut inputs.constraints, leaf_key, constraint.clone());
    set_edge_branch_dist(&graph, &mut inputs, leaf_key, 3.0);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let actual = node_time_distribution(&backward, leaf_key).expect("leaf A should have a time distribution");
    assert_eq!(&constraint, actual.as_ref());

    Ok(())
  }

  #[test]
  fn test_backward_pass_applies_internal_node_date_constraint() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal I not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_date_constraint(
      &mut inputs.constraints,
      leaf_a_key,
      Distribution::range((2014.0, 2016.0), 0.0),
    );
    set_date_constraint(
      &mut inputs.constraints,
      leaf_b_key,
      Distribution::range((2013.0, 2015.0), 0.0),
    );
    set_date_constraint(
      &mut inputs.constraints,
      internal_key,
      Distribution::range((2012.0, 2014.0), 0.0),
    );
    set_edge_branch_dist(&graph, &mut inputs, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_b_key, 2.0);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let actual =
      node_time_distribution(&backward, internal_key).expect("internal node should have a time distribution");
    let expected = Distribution::range((2012.0, 2013.0), 0.0);
    assert_eq!(&expected, actual.as_ref());

    Ok(())
  }

  #[test]
  fn test_backward_pass_sets_edge_messages() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:2.5)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_key, 2013.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_key, 2.5);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    for edge in graph.get_edges() {
      let edge_read = edge;
      if edge_read.target() == leaf_key {
        let msg = backward.messages[&edge_read.key()]
          .as_ref()
          .expect("edge should have msg_to_parent after backward pass");
        let msg_time = msg.likely_time()?.expect("message should have likely_time");
        pretty_assert_ulps_eq!(msg_time, 2010.5, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_skips_bad_branch_children() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_a_key, 2015.0);
    set_leaf_time(&mut inputs, leaf_b_key, 2014.0);

    inputs.bad_branches.insert(leaf_b_key, true);

    set_edge_branch_dist(&graph, &mut inputs, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_b_key, 2.0);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let time_dist =
      node_time_distribution(&backward, internal_key).expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()?
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_backward_pass_bad_branch_equivalent_to_removal() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0)I:1.0)root;")?;
    let ref_names = nwk_parsed.names();
    let ref_graph = nwk_parsed.graph;
    let ref_a_key = find_node_key_by_name(&ref_graph, &ref_names, "A").expect("leaf A not found");
    let mut ref_inputs = BackwardInputs::new(&ref_graph);
    set_leaf_time(&mut ref_inputs, ref_a_key, 2015.0);
    set_edge_branch_dist(&ref_graph, &mut ref_inputs, ref_a_key, 3.0);
    let ref_backward = run_backward_pass(&ref_graph, &ref_inputs, None)?;

    let ref_internal_key = find_node_key_by_name(&ref_graph, &ref_names, "I").expect("internal I not found");
    let ref_time = node_time_distribution(&ref_backward, ref_internal_key)
      .expect("should have time dist")
      .likely_time()?
      .expect("should have likely_time");

    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let test_names = nwk_parsed.names();
    let test_graph = nwk_parsed.graph;
    let test_a_key = find_node_key_by_name(&test_graph, &test_names, "A").expect("leaf A not found");
    let test_b_key = find_node_key_by_name(&test_graph, &test_names, "B").expect("leaf B not found");
    let mut test_inputs = BackwardInputs::new(&test_graph);
    set_leaf_time(&mut test_inputs, test_a_key, 2015.0);
    set_leaf_time(&mut test_inputs, test_b_key, 2014.0);
    set_edge_branch_dist(&test_graph, &mut test_inputs, test_a_key, 3.0);
    set_edge_branch_dist(&test_graph, &mut test_inputs, test_b_key, 2.0);

    test_inputs.bad_branches.insert(test_b_key, true);

    let test_backward = run_backward_pass(&test_graph, &test_inputs, None)?;

    let test_internal_key = find_node_key_by_name(&test_graph, &test_names, "I").expect("internal I not found");
    let test_time = node_time_distribution(&test_backward, test_internal_key)
      .expect("should have time dist")
      .likely_time()?
      .expect("should have likely_time");

    pretty_assert_ulps_eq!(ref_time, test_time, max_ulps = 4);

    Ok(())
  }

  #[test]
  #[ignore = "fold now receives mass-windowed messages; Gaussian-product oracle no longer holds"]
  fn test_backward_pass_sums_function_children_to_gaussian_product() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let a = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let b = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let c = find_node_key_by_name(&graph, &names, "C").expect("leaf C not found");

    let x = Array1::linspace(2000.0, 2010.0, 2001);
    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_function(&mut inputs, a, &x, gaussian_neglog(&x, 2002.0, 1.0))?;
    set_leaf_function(&mut inputs, b, &x, gaussian_neglog(&x, 2008.0, 1.0))?;
    set_leaf_function(&mut inputs, c, &x, gaussian_neglog(&x, 2005.0, 2.0))?;
    set_edge_branch_dist(&graph, &mut inputs, a, 0.0);
    set_edge_branch_dist(&graph, &mut inputs, b, 0.0);
    set_edge_branch_dist(&graph, &mut inputs, c, 0.0);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let internal = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let dist = node_time_distribution(&backward, internal).expect("internal node should have a time distribution");

    let grid = dist.t();
    let spacing = grid[1] - grid[0];
    let peak = dist.likely_time()?.expect("distribution should have a likely_time");
    assert_abs_diff_eq!(peak, 2005.0, epsilon = spacing);

    for t in [2004.0_f64, 2005.0, 2006.0] {
      let expected = 2.0 * (t - 2005.0).powi(2);
      assert_abs_diff_eq!(dist.eval(t)?, expected, epsilon = 1e-4);
    }

    Ok(())
  }

  #[test]
  fn test_backward_pass_fan_out_result_independent_of_child_order() -> Result<(), Report> {
    let x = Array1::linspace(2000.0, 2010.0, 11);
    let ya = gaussian_neglog(&x, 2002.0, 1.0);
    let yb = gaussian_neglog(&x, 2008.0, 1.0);
    let yc = gaussian_neglog(&x, 2005.0, 2.0);

    let fold_in_order = |newick: &str| -> Result<(Array1<f64>, f64), Report> {
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let mut inputs = BackwardInputs::new(&graph);
      for (name, y) in [("A", &ya), ("B", &yb), ("C", &yc)] {
        let key = find_node_key_by_name(&graph, &names, name).expect("leaf not found");
        set_leaf_function(&mut inputs, key, &x, y.clone())?;
        set_edge_branch_dist(&graph, &mut inputs, key, 0.0);
      }
      let backward = run_backward_pass(&graph, &inputs, None)?;
      let internal = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
      let dist = node_time_distribution(&backward, internal).expect("internal node should have a time distribution");
      Ok((
        dist.y()?,
        dist.likely_time()?.expect("distribution should have a likely_time"),
      ))
    };

    let (y_abc, peak_abc) = fold_in_order("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let (y_cab, peak_cab) = fold_in_order("((C:0.0,A:0.0,B:0.0)I:1.0)root;")?;

    pretty_assert_ulps_eq!(y_abc, y_cab, max_ulps = 8);
    pretty_assert_ulps_eq!(peak_abc, peak_cab, max_ulps = 4);

    Ok(())
  }

  #[test]
  #[ignore = "fold now receives mass-windowed messages; precision-weighted-mean oracle no longer holds"]
  fn test_backward_pass_function_children_peak_at_precision_weighted_mean() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let a = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let b = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let c = find_node_key_by_name(&graph, &names, "C").expect("leaf C not found");

    let x = Array1::linspace(2000.0, 2010.0, 11);
    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_function(&mut inputs, a, &x, gaussian_neglog(&x, 2002.0, 1.0))?;
    set_leaf_function(&mut inputs, b, &x, gaussian_neglog(&x, 2008.0, 1.0))?;
    set_leaf_function(&mut inputs, c, &x, gaussian_neglog(&x, 2005.0, 2.0))?;
    set_edge_branch_dist(&graph, &mut inputs, a, 0.0);
    set_edge_branch_dist(&graph, &mut inputs, b, 0.0);
    set_edge_branch_dist(&graph, &mut inputs, c, 0.0);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let internal = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let dist = node_time_distribution(&backward, internal).expect("internal node should have a time distribution");
    let likely_time = dist.likely_time()?.expect("distribution should have a likely_time");

    let grid = dist.t();
    let spacing = grid[1] - grid[0];
    assert_abs_diff_eq!(likely_time, 2005.0, epsilon = spacing);

    Ok(())
  }

  #[test]
  fn test_backward_pass_bad_branch_sends_no_message() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_a_key, 2015.0);
    set_leaf_time(&mut inputs, leaf_b_key, 2014.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_b_key, 2.0);
    inputs.bad_branches.insert(leaf_b_key, true);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    assert_eq!(None, backward.messages[&parent_edge_key(&graph, leaf_b_key)]);
    assert!(backward.messages[&parent_edge_key(&graph, leaf_a_key)].is_some());
    Ok(())
  }

  #[test]
  fn test_backward_pass_node_without_evidence_has_no_subtree_distribution() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:3.0,B:2.0)I:1.0,C:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let leaf_c_key = find_node_key_by_name(&graph, &names, "C").expect("leaf C not found");
    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal I not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_c_key, 2015.0);
    for key in [leaf_a_key, leaf_b_key, internal_key, leaf_c_key] {
      set_edge_branch_dist(&graph, &mut inputs, key, 1.0);
    }
    inputs.bad_branches.insert(leaf_a_key, true);
    inputs.bad_branches.insert(leaf_b_key, true);
    inputs.bad_branches.insert(internal_key, true);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    assert_eq!(None, node_time_distribution(&backward, internal_key));
    assert_eq!(None, backward.messages[&parent_edge_key(&graph, internal_key)]);
    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;

    pub(super) struct BackwardInputs {
      pub constraints: DateConstraints,
      pub bad_branches: BTreeMap<GraphNodeKey, bool>,
      pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
    }

    impl BackwardInputs {
      pub(super) fn new(graph: &Graph) -> Self {
        let bad_branches = graph.get_nodes().map(|node| (node.key(), false)).collect();
        let branches = graph
          .get_edges()
          .map(|edge| {
            let branch = BranchLikelihood {
              distribution: None,
              time_length: None,
            };
            (edge.key(), branch)
          })
          .collect();
        Self {
          constraints: DateConstraints::default(),
          bad_branches,
          branches,
        }
      }
    }

    pub(super) fn run_backward_pass(
      graph: &Graph,
      inputs: &BackwardInputs,
      coalescent_model: Option<&CoalescentModel>,
    ) -> Result<TimeBackward, Report> {
      propagate_distributions_backward(
        graph,
        &inputs.constraints,
        coalescent_model,
        &inputs.bad_branches,
        &inputs.branches,
      )
    }

    pub(super) fn node_time_distribution(
      backward: &TimeBackward,
      key: GraphNodeKey,
    ) -> Option<Arc<Distribution<NegLog>>> {
      backward.subtree[&key].clone()
    }

    pub(super) fn parent_edge_key(graph: &Graph, target_key: GraphNodeKey) -> GraphEdgeKey {
      graph
        .get_edges()
        .find(|edge| edge.target() == target_key)
        .expect("node must have a parent edge")
        .key()
    }

    pub(super) fn set_date_constraint(
      constraints: &mut DateConstraints,
      key: GraphNodeKey,
      dist: Distribution<NegLog>,
    ) {
      constraints.date_constraints.insert(key, Some(Arc::new(dist)));
    }

    pub(super) fn set_leaf_time(inputs: &mut BackwardInputs, key: GraphNodeKey, time: f64) {
      set_date_constraint(&mut inputs.constraints, key, Distribution::point(time, 0.0));
    }

    pub(super) fn set_edge_branch_dist(graph: &Graph, inputs: &mut BackwardInputs, target_key: GraphNodeKey, bl: f64) {
      let branch = BranchLikelihood {
        distribution: Some(Arc::new(Distribution::point(bl, 0.0))),
        time_length: Some(bl),
      };
      inputs.branches.insert(parent_edge_key(graph, target_key), branch);
    }

    pub(super) fn gaussian_neglog(x: &Array1<f64>, mean: f64, precision: f64) -> Array1<f64> {
      x.mapv(|t| 0.5 * precision * (t - mean).powi(2))
    }

    pub(super) fn set_leaf_function(
      inputs: &mut BackwardInputs,
      key: GraphNodeKey,
      x: &Array1<f64>,
      y: Array1<f64>,
    ) -> Result<(), Report> {
      let dist = Distribution::function(x.clone(), y)?;
      set_date_constraint(&mut inputs.constraints, key, dist);
      Ok(())
    }

    pub(super) fn coalescent_model(tc: f64) -> Result<CoalescentModel, Report> {
      let lineage_counts = PiecewiseConstantFn::new(array![1900.0, 2100.0], array![1.0, 2.0, 0.0]);
      CoalescentModel::new(&lineage_counts, &Distribution::constant(tc))
    }
  }

  use helpers::*;
}
