#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, unknown_branches};
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::result::{BranchLikelihood, TimeBackward};
  use eyre::Report;
  use ndarray::Array1;
  use ndarray::array;
  use ordered_float::OrderedFloat;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;
  use treetime_utils::pretty_assert_abs_diff_eq;

  const LAPLACE_WEIGHTED_MEDIAN: f64 = 2005.0;

  const LAPLACE_PROBE_NEAR: f64 = 2006.0;

  const LAPLACE_PROBE_FAR: f64 = 2007.0;

  const BRANCH_GRID_POINTS: usize = 201;

  const BRANCH_PRECISION: f64 = 4.0;

  #[test]
  fn test_backward_pass_computes_internal_node_time() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:2.5)I:1.0)root;".as_slice())?;
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
    let nwk_parsed = nwk_read(b"((A:2.5)I:1.0)root;".as_slice())?;
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
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
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
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
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

    let backward = run_backward_pass(&graph, &inputs, Some(&coalescent_model))?;

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
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
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
  fn test_backward_pass_coalescent_prior_adds_the_internal_merger_term_to_the_subtree() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let internal_key = find_node_key_by_name(&graph, &names, "I").expect("internal I not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_a_key, 2015.0);
    set_leaf_time(&mut inputs, leaf_b_key, 2014.0);
    set_edge_branch_function(&graph, &mut inputs, leaf_a_key, 3.0)?;
    set_edge_branch_function(&graph, &mut inputs, leaf_b_key, 2.0)?;
    let model = coalescent_model(0.5)?;

    let without = run_backward_pass(&graph, &inputs, None)?;
    let with = run_backward_pass(&graph, &inputs, Some(&model))?;

    let without = node_time_distribution(&without, internal_key).expect("I must have a subtree distribution");
    let with = node_time_distribution(&with, internal_key).expect("I must have a subtree distribution");
    assert_eq!(without.t(), with.t());
    let grid = without.t();
    let reference = grid[0];
    let n_children = 2;
    for &time in grid.iter().skip(1) {
      let added = (with.eval(time)? - with.eval(reference)?) - (without.eval(time)? - without.eval(reference)?);
      let expected =
        model.internal_contribution(time, n_children)? - model.internal_contribution(reference, n_children)?;
      pretty_assert_abs_diff_eq!(added, expected, epsilon = 1e-9);
    }
    Ok(())
  }

  #[test]
  fn test_backward_pass_leaf_subtree_distribution_is_its_date_constraint() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0)I:1.0)root;".as_slice())?;
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
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
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
    let nwk_parsed = nwk_read(b"((A:2.5)I:1.0)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;

    let leaf_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_key, 2013.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_key, 2.5);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    let msg = backward.messages[&parent_edge_key(&graph, leaf_key)]
      .as_ref()
      .expect("edge should have msg_to_parent after backward pass");
    let msg_time = msg.likely_time()?.expect("message should have likely_time");
    pretty_assert_ulps_eq!(msg_time, 2010.5, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_backward_pass_skips_bad_branch_children() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
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
    let nwk_parsed = nwk_read(b"((A:3.0)I:1.0)root;".as_slice())?;
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

    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
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
  fn test_backward_pass_sums_function_children_neglogs() -> Result<(), Report> {
    let (dist, _) = fold_laplace_children()?;

    let slope_right_of_the_median = dist.eval(LAPLACE_PROBE_FAR)? - dist.eval(LAPLACE_PROBE_NEAR)?;

    pretty_assert_abs_diff_eq!(slope_right_of_the_median, 2.0, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_backward_pass_function_children_peak_at_the_weighted_median() -> Result<(), Report> {
    let (dist, peak) = fold_laplace_children()?;

    let grid = dist.t();
    let nearest_to_the_median = grid
      .iter()
      .copied()
      .min_by_key(|time| OrderedFloat((time - LAPLACE_WEIGHTED_MEDIAN).abs()))
      .expect("the folded distribution must have a grid");
    pretty_assert_ulps_eq!(nearest_to_the_median, peak, max_ulps = 0);
    Ok(())
  }

  #[test]
  fn test_backward_pass_fan_out_result_independent_of_child_order() -> Result<(), Report> {
    let x = Array1::linspace(2000.0, 2010.0, 11);
    let ya = gaussian_neglog(&x, 2002.0, 1.0);
    let yb = gaussian_neglog(&x, 2008.0, 1.0);
    let yc = gaussian_neglog(&x, 2005.0, 2.0);

    let fold_in_order = |newick: &str| -> Result<(Array1<f64>, f64), Report> {
      let nwk_parsed = nwk_read(newick.as_bytes())?;
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
  fn test_backward_pass_bad_branch_sends_no_message() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
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
  fn test_backward_pass_sends_no_message_without_a_subtree_or_a_branch_likelihood() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0,C:1.0)I:1.0)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, &names, "B").expect("leaf B not found");
    let leaf_c_key = find_node_key_by_name(&graph, &names, "C").expect("leaf C not found");

    let mut inputs = BackwardInputs::new(&graph);
    set_leaf_time(&mut inputs, leaf_a_key, 2015.0);
    set_leaf_time(&mut inputs, leaf_b_key, 2014.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, &mut inputs, leaf_c_key, 1.0);

    let backward = run_backward_pass(&graph, &inputs, None)?;

    assert!(backward.messages[&parent_edge_key(&graph, leaf_a_key)].is_some());
    assert_eq!(None, backward.messages[&parent_edge_key(&graph, leaf_b_key)]);
    assert_eq!(None, node_time_distribution(&backward, leaf_c_key));
    assert_eq!(None, backward.messages[&parent_edge_key(&graph, leaf_c_key)]);
    Ok(())
  }

  #[test]
  fn test_backward_pass_missing_bad_branch_flag_is_an_internal_error() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0)root;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let leaf_a_key = find_node_key_by_name(&graph, &names, "A").expect("leaf A not found");

    let mut inputs = BackwardInputs::new(&graph);
    inputs.bad_branches.remove(&leaf_a_key);

    assert_error!(
      run_backward_pass(&graph, &inputs, None),
      format!(
        "Bad-branch flags are missing node {leaf_a_key}. This is an internal error. Please report it to developers."
      )
    );
    Ok(())
  }

  #[test]
  fn test_backward_pass_node_without_evidence_has_no_subtree_distribution() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:3.0,B:2.0)I:1.0,C:1.0)root;".as_slice())?;
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
        Self {
          constraints: DateConstraints::default(),
          bad_branches,
          branches: unknown_branches(graph),
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
        MaxGridPoints::default(),
      )
    }

    pub(super) fn node_time_distribution(
      backward: &TimeBackward,
      key: GraphNodeKey,
    ) -> Option<Arc<Distribution<NegLog>>> {
      backward.subtree[&key].clone()
    }

    pub(super) fn set_date_constraint(
      constraints: &mut DateConstraints,
      key: GraphNodeKey,
      dist: Distribution<NegLog>,
    ) {
      constraints.by_node.insert(key, Some(Arc::new(dist)));
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

    pub(super) fn set_edge_branch_function(
      graph: &Graph,
      inputs: &mut BackwardInputs,
      target_key: GraphNodeKey,
      mean: f64,
    ) -> Result<(), Report> {
      let t = Array1::linspace(0.0, 2.0 * mean, BRANCH_GRID_POINTS);
      let y = gaussian_neglog(&t, mean, BRANCH_PRECISION);
      let branch = BranchLikelihood {
        distribution: Some(Arc::new(Distribution::function(t, y)?)),
        time_length: Some(mean),
      };
      inputs.branches.insert(parent_edge_key(graph, target_key), branch);
      Ok(())
    }

    pub(super) fn laplace_neglog(x: &Array1<f64>, median: f64, weight: f64) -> Array1<f64> {
      x.mapv(|t| weight * (t - median).abs())
    }

    pub(super) fn fold_laplace_children() -> Result<(Arc<Distribution<NegLog>>, f64), Report> {
      let nwk_parsed = nwk_read(b"((A:0.0,B:0.0,C:0.0)I:1.0)root;".as_slice())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let x = Array1::linspace(2000.0, 2010.0, 2001);
      let mut inputs = BackwardInputs::new(&graph);
      for (name, median, weight) in [("A", 2002.0, 1.0), ("B", 2008.0, 1.0), ("C", 2005.0, 2.0)] {
        let key = find_node_key_by_name(&graph, &names, name).expect("leaf not found");
        set_leaf_function(&mut inputs, key, &x, laplace_neglog(&x, median, weight))?;
        set_edge_branch_dist(&graph, &mut inputs, key, 0.0);
      }

      let backward = run_backward_pass(&graph, &inputs, None)?;

      let internal = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
      let dist = node_time_distribution(&backward, internal).expect("internal node should have a time distribution");
      let peak = dist.likely_time()?.expect("distribution should have a likely_time");
      Ok((dist, peak))
    }

    pub(super) fn coalescent_model(tc: f64) -> Result<CoalescentModel, Report> {
      let lineage_counts = PiecewiseConstantFn::new(array![1900.0, 2100.0], array![1.0, 2.0, 0.0]);
      CoalescentModel::new(&lineage_counts, &Distribution::constant(tc))
    }
  }

  use helpers::*;
}
