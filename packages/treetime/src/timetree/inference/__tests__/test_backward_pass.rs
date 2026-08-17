#[cfg(test)]
mod tests {
  use crate::coalescent::coalescent::CoalescentModel;
  use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use ndarray::array;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::edge::BranchDistribution;
  use treetime_graph::node::{GraphNodeKey, TimeConstraint};
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::nwk::nwk_read_str;

  /// Test backward pass on a simple 2-leaf tree:
  /// ((A:2.5)I:1.0)root;
  /// A has time 2013.0, I should get 2013.0 - 2.5 = 2010.5
  #[test]
  fn test_backward_pass_computes_internal_node_time() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:2.5)I:1.0)root;")?;

    // Set leaf A's time distribution to point at 2013.0
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    {
      let leaf_node = graph.get_node(leaf_key).expect("leaf A exists");
      let mut payload = leaf_node.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2013.0, 0.0)));
    }

    // Set branch length distribution on edge from I to A (branch length 2.5 years)
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      if edge_read.target() == leaf_key {
        let mut payload = edge_read.payload().write_arc();
        payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(2.5, 0.0))));
      }
    }

    // Run backward pass
    propagate_distributions_backward(&graph, None)?;

    // Check internal node I has time distribution centered at 2013.0 - 2.5 = 2010.5
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let internal_node = graph.get_node(internal_key).expect("internal node I exists");
    let payload = internal_node.read_arc().payload().read_arc();

    let time_dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have time distribution after backward pass");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2010.5, max_ulps = 4);
    assert!(likely_time < 2013.0, "Parent should be older than child");

    Ok(())
  }

  /// Test backward pass with two children: messages are multiplied (intersection).
  /// Tree: ((A:3.0,B:2.0)I:1.0)root;
  /// A at 2015.0, B at 2014.0
  /// I gets messages: from A -> 2015-3=2012, from B -> 2014-2=2012
  /// Both agree, so I should be at 2012.0
  #[test]
  fn test_backward_pass_multiplies_child_messages() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;

    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");

    // Set time distributions on leaves
    {
      let node_a = graph.get_node(leaf_a_key).expect("leaf A exists");
      let mut payload = node_a.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2015.0, 0.0)));
    }
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      let mut payload = node_b.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2014.0, 0.0)));
    }

    // Set branch length distributions
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      let target = edge_read.target();
      let branch_length = if target == leaf_a_key {
        3.0
      } else if target == leaf_b_key {
        2.0
      } else {
        continue;
      };
      let mut payload = edge_read.payload().write_arc();
      payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(branch_length, 0.0))));
    }

    propagate_distributions_backward(&graph, None)?;

    // Internal node I should have time at 2012.0 (both children agree)
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let internal_node = graph.get_node(internal_key).expect("internal node I exists");
    let payload = internal_node.read_arc().payload().read_arc();

    let time_dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  /// Test that backward pass preserves leaf time_distribution when coalescent is provided.
  /// This is a regression test for a bug where coalescent contributions overwrote leaf dates,
  /// causing subsequent clock regression to fail with "No variation in sampling dates".
  #[test]
  fn test_backward_pass_preserves_leaf_time_distribution_with_coalescent() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;

    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");

    let date_a = 2015.0;
    let date_b = 2014.0;

    // Set time distributions on leaves (date constraints)
    {
      let node_a = graph.get_node(leaf_a_key).expect("leaf A exists");
      let mut payload = node_a.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(date_a, 0.0)));
    }
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      let mut payload = node_b.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(date_b, 0.0)));
    }

    // Set branch length distributions
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      let target = edge_read.target();
      let branch_length = if target == leaf_a_key {
        3.0
      } else if target == leaf_b_key {
        2.0
      } else {
        continue;
      };
      let mut payload = edge_read.payload().write_arc();
      payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(branch_length, 0.0))));
    }

    let coalescent_model = coalescent_model(0.01)?;

    propagate_distributions_backward(&graph, Some(&coalescent_model))?;
    propagate_distributions_backward(&graph, Some(&coalescent_model))?;

    // Verify leaf A still has its original date
    {
      let node_a = graph.get_node(leaf_a_key).expect("leaf A exists");
      let payload = node_a.read_arc().payload().read_arc();
      let time_dist = payload
        .time_distribution()
        .as_ref()
        .expect("leaf A should have time distribution");
      let expected = Distribution::point(date_a, 0.0);
      assert_eq!(&expected, time_dist.as_ref());
    }

    // Verify leaf B still has its original date
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      let payload = node_b.read_arc().payload().read_arc();
      let time_dist = payload
        .time_distribution()
        .as_ref()
        .expect("leaf B should have time distribution");
      let expected = Distribution::point(date_b, 0.0);
      assert_eq!(&expected, time_dist.as_ref());
    }

    Ok(())
  }

  /// V0 normalizes negative-log likelihoods relative to their peak before
  /// exponentiation (`treetime/distribution.py:210-213`). A strong coalescent
  /// factor must therefore preserve a point-supported internal date.
  #[test]
  fn test_backward_pass_preserves_internal_time_with_strong_coalescent() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;
    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal I not found");

    set_leaf_time(&graph, leaf_a_key, 2015.0);
    set_leaf_time(&graph, leaf_b_key, 2014.0);
    set_edge_branch_dist(&graph, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, leaf_b_key, 2.0);

    let coalescent_model = coalescent_model(1e-6)?;

    propagate_distributions_backward(&graph, Some(&coalescent_model))?;

    let internal = graph.get_node(internal_key).expect("internal I exists");
    let payload = internal.read_arc().payload().read_arc();
    let actual = payload
      .time_distribution()
      .as_ref()
      .and_then(|distribution| distribution.likely_time());
    let expected = Some(2012.0);
    assert_eq!(expected, actual);

    Ok(())
  }

  /// The date given as input is fixed for the whole run, so each backward pass lifts it back into
  /// the time distribution. The forward pass refines the time distribution of a leaf whose date is
  /// uncertain in place; without the lift that refined distribution would be sent up to the parent
  /// on the next round as if it were an independent observation.
  #[test]
  fn test_backward_pass_restores_leaf_time_distribution_from_the_date_constraint() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0)I:1.0)root;")?;
    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    let constraint = Distribution::range((2014.0, 2015.0), 0.0);
    {
      let node = graph.get_node(leaf_key).expect("leaf A exists");
      let mut payload = node.read_arc().payload().write_arc();
      payload.date_constraint = Some(Arc::new(constraint.clone()));
      // What a forward pass leaves behind: the range narrowed down by the rest of the tree.
      payload.time_distribution = Some(Arc::new(Distribution::point(2014.2, 0.0)));
    }
    set_edge_branch_dist(&graph, leaf_key, 3.0);

    propagate_distributions_backward(&graph, None)?;

    let payload = graph
      .get_node(leaf_key)
      .expect("leaf A exists")
      .read_arc()
      .payload()
      .read_arc();
    let actual = payload
      .time_distribution()
      .as_ref()
      .expect("leaf A should have a time distribution");
    assert_eq!(&constraint, actual.as_ref());

    Ok(())
  }

  /// A date given for an internal node constrains it the same way a leaf date does: it multiplies
  /// what the children have to say rather than being replaced by it. Here both children put the
  /// node somewhere in [2011, 2013] and its own date narrows that to [2012, 2013].
  #[test]
  fn test_backward_pass_applies_internal_node_date_constraint() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;
    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal I not found");

    set_date_constraint(&graph, leaf_a_key, Distribution::range((2014.0, 2016.0), 0.0));
    set_date_constraint(&graph, leaf_b_key, Distribution::range((2013.0, 2015.0), 0.0));
    set_date_constraint(&graph, internal_key, Distribution::range((2012.0, 2014.0), 0.0));
    set_edge_branch_dist(&graph, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, leaf_b_key, 2.0);

    propagate_distributions_backward(&graph, None)?;

    let payload = graph
      .get_node(internal_key)
      .expect("internal I exists")
      .read_arc()
      .payload()
      .read_arc();
    let actual = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have a time distribution");
    let expected = Distribution::range((2012.0, 2013.0), 0.0);
    assert_eq!(&expected, actual.as_ref());

    Ok(())
  }

  /// Test that backward pass stores msg_to_parent on edges.
  #[test]
  fn test_backward_pass_sets_edge_messages() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:2.5)I:1.0)root;")?;

    let leaf_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");

    // Set leaf time distribution
    {
      let leaf_node = graph.get_node(leaf_key).expect("leaf A exists");
      let mut payload = leaf_node.read_arc().payload().write_arc();
      payload.time_distribution = Some(Arc::new(Distribution::point(2013.0, 0.0)));
    }

    // Set branch length distribution
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      if edge_read.target() == leaf_key {
        let mut payload = edge_read.payload().write_arc();
        payload.set_branch_length_distribution(Some(Arc::new(Distribution::point(2.5, 0.0))));
      }
    }

    propagate_distributions_backward(&graph, None)?;

    // Check edge from I to A has msg_to_parent set
    for edge in graph.get_edges() {
      let edge_read = edge.read_arc();
      if edge_read.target() == leaf_key {
        let payload = edge_read.payload().read_arc();
        let msg = payload
          .msg_to_parent()
          .as_ref()
          .expect("edge should have msg_to_parent after backward pass");
        let msg_time = msg.likely_time().expect("message should have likely_time");
        // Message should be parent time: 2013.0 - 2.5 = 2010.5
        pretty_assert_ulps_eq!(msg_time, 2010.5, max_ulps = 4);
      }
    }

    Ok(())
  }

  /// Test backward pass skips children with bad_branch=true.
  /// Tree: ((A:3.0,B:2.0)I:1.0)root;
  /// A at 2015.0, B at 2014.0. Mark B as bad_branch.
  /// I should get time from A only: 2015.0 - 3.0 = 2012.0
  #[test]
  fn test_backward_pass_skips_bad_branch_children() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;

    let leaf_a_key = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let leaf_b_key = find_node_key_by_name(&graph, "B").expect("leaf B not found");

    // Set time distributions on leaves
    set_leaf_time(&graph, leaf_a_key, 2015.0);
    set_leaf_time(&graph, leaf_b_key, 2014.0);

    // Mark B as bad_branch
    {
      let node_b = graph.get_node(leaf_b_key).expect("leaf B exists");
      node_b.read_arc().payload().write_arc().bad_branch = true;
    }

    // Set branch length distributions on both edges
    set_edge_branch_dist(&graph, leaf_a_key, 3.0);
    set_edge_branch_dist(&graph, leaf_b_key, 2.0);

    propagate_distributions_backward(&graph, None)?;

    // I should get time only from A: 2015.0 - 3.0 = 2012.0
    let internal_key = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let internal_node = graph.get_node(internal_key).expect("internal node I exists");
    let payload = internal_node.read_arc().payload().read_arc();

    let time_dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have time distribution");
    let likely_time = time_dist
      .likely_time()
      .expect("time distribution should have likely_time");

    pretty_assert_ulps_eq!(likely_time, 2012.0, max_ulps = 4);

    Ok(())
  }

  /// Test that marking a leaf as bad_branch produces identical parent time
  /// to a tree where that leaf is absent.
  #[test]
  fn test_backward_pass_bad_branch_equivalent_to_removal() -> Result<(), Report> {
    // Reference tree: only A
    let ref_graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0)I:1.0)root;")?;
    let ref_a_key = find_node_key_by_name(&ref_graph, "A").expect("leaf A not found");
    set_leaf_time(&ref_graph, ref_a_key, 2015.0);
    set_edge_branch_dist(&ref_graph, ref_a_key, 3.0);
    propagate_distributions_backward(&ref_graph, None)?;

    let ref_internal_key = find_node_key_by_name(&ref_graph, "I").expect("internal I not found");
    let ref_time = ref_graph
      .get_node(ref_internal_key)
      .expect("I exists")
      .read_arc()
      .payload()
      .read_arc()
      .time_distribution()
      .as_ref()
      .expect("should have time dist")
      .likely_time()
      .expect("should have likely_time");

    // Test tree: A + B(bad)
    let test_graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:3.0,B:2.0)I:1.0)root;")?;
    let test_a_key = find_node_key_by_name(&test_graph, "A").expect("leaf A not found");
    let test_b_key = find_node_key_by_name(&test_graph, "B").expect("leaf B not found");
    set_leaf_time(&test_graph, test_a_key, 2015.0);
    set_leaf_time(&test_graph, test_b_key, 2014.0);
    set_edge_branch_dist(&test_graph, test_a_key, 3.0);
    set_edge_branch_dist(&test_graph, test_b_key, 2.0);

    // Mark B as bad
    test_graph
      .get_node(test_b_key)
      .expect("B exists")
      .read_arc()
      .payload()
      .write_arc()
      .bad_branch = true;

    propagate_distributions_backward(&test_graph, None)?;

    let test_internal_key = find_node_key_by_name(&test_graph, "I").expect("internal I not found");
    let test_time = test_graph
      .get_node(test_internal_key)
      .expect("I exists")
      .read_arc()
      .payload()
      .read_arc()
      .time_distribution()
      .as_ref()
      .expect("should have time dist")
      .likely_time()
      .expect("should have likely_time");

    pretty_assert_ulps_eq!(ref_time, test_time, max_ulps = 4);

    Ok(())
  }

  /// Three children whose backward messages are gridded `Function` distributions on a shared grid.
  /// A zero-length branch makes each message equal to its child's distribution, so the fold reduces
  /// to a pure sum of neg-log ordinates. The product of the three Gaussians is itself Gaussian, with
  /// precision `1 + 1 + 2 = 4` and mean `(2002 + 2008 + 2*2005)/4 = 2005`, so its peak-normalized
  /// neg-log is the parabola `2*(t - 2005)^2` -- an independent analytic oracle for the fold. Mass
  /// re-windowing then regrids that posterior onto its mass domain, so the check is against the
  /// analytic function rather than the child grid.
  #[test]
  fn test_backward_pass_sums_function_children_to_gaussian_product() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let a = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let b = find_node_key_by_name(&graph, "B").expect("leaf B not found");
    let c = find_node_key_by_name(&graph, "C").expect("leaf C not found");

    // A fine child grid keeps the re-window's linear resampling faithful to the analytic parabola.
    let x = Array1::linspace(2000.0, 2010.0, 2001);
    set_leaf_function(&graph, a, &x, gaussian_neglog(&x, 2002.0, 1.0))?;
    set_leaf_function(&graph, b, &x, gaussian_neglog(&x, 2008.0, 1.0))?;
    set_leaf_function(&graph, c, &x, gaussian_neglog(&x, 2005.0, 2.0))?;
    set_edge_branch_dist(&graph, a, 0.0);
    set_edge_branch_dist(&graph, b, 0.0);
    set_edge_branch_dist(&graph, c, 0.0);

    propagate_distributions_backward(&graph, None)?;

    let internal = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let node = graph.get_node(internal).expect("internal I exists");
    let payload = node.read_arc().payload().read_arc();
    let dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have a time distribution");

    // Peak sits at the precision-weighted mean, within one grid spacing after re-windowing.
    let grid = dist.t();
    let spacing = grid[1] - grid[0];
    let peak = dist.likely_time().expect("distribution should have a likely_time");
    assert_abs_diff_eq!(peak, 2005.0, epsilon = spacing);

    // Peak-normalized neg-log matches the analytic Gaussian-product parabola at interior points. The
    // tolerance covers the O(dx^2) linear-resampling error and the peak-normalization offset (the
    // stored minimum sits within dx/2 of the true mode 2005); both are ~1e-5 at this resolution.
    for t in [2004.0_f64, 2005.0, 2006.0] {
      let expected = 2.0 * (t - 2005.0).powi(2);
      assert_abs_diff_eq!(dist.eval(t)?, expected, epsilon = 1e-4);
    }

    Ok(())
  }

  /// A fan-out node's folded distribution must not depend on the order its children are visited.
  /// The same three Function messages folded from two different child orderings must agree; the old
  /// per-child multiply-and-resample fold did not, because each step resampled the accumulator.
  #[test]
  fn test_backward_pass_fan_out_result_independent_of_child_order() -> Result<(), Report> {
    let x = Array1::linspace(2000.0, 2010.0, 11);
    let ya = gaussian_neglog(&x, 2002.0, 1.0);
    let yb = gaussian_neglog(&x, 2008.0, 1.0);
    let yc = gaussian_neglog(&x, 2005.0, 2.0);

    let fold_in_order = |newick: &str| -> Result<(Array1<f64>, f64), Report> {
      let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>(newick)?;
      for (name, y) in [("A", &ya), ("B", &yb), ("C", &yc)] {
        let key = find_node_key_by_name(&graph, name).expect("leaf not found");
        set_leaf_function(&graph, key, &x, y.clone())?;
        set_edge_branch_dist(&graph, key, 0.0);
      }
      propagate_distributions_backward(&graph, None)?;
      let internal = find_node_key_by_name(&graph, "I").expect("internal node I not found");
      let node = graph.get_node(internal).expect("internal I exists");
      let payload = node.read_arc().payload().read_arc();
      let dist = payload
        .time_distribution()
        .as_ref()
        .expect("internal node should have a time distribution");
      Ok((
        dist.y(),
        dist.likely_time().expect("distribution should have a likely_time"),
      ))
    };

    let (y_abc, peak_abc) = fold_in_order("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let (y_cab, peak_cab) = fold_in_order("((C:0.0,A:0.0,B:0.0)I:1.0)root;")?;

    pretty_assert_ulps_eq!(y_abc, y_cab, max_ulps = 8);
    pretty_assert_ulps_eq!(peak_abc, peak_cab, max_ulps = 4);

    Ok(())
  }

  /// The product of Gaussian likelihoods is Gaussian with the precision-weighted mean of the
  /// operands. In neg-log space each message is a parabola, and the folded peak must sit at that
  /// weighted mean -- an independent analytic oracle for the summed fold.
  #[test]
  fn test_backward_pass_function_children_peak_at_precision_weighted_mean() -> Result<(), Report> {
    let graph = nwk_read_str::<NodeTimetree, EdgeTimetree, ()>("((A:0.0,B:0.0,C:0.0)I:1.0)root;")?;
    let a = find_node_key_by_name(&graph, "A").expect("leaf A not found");
    let b = find_node_key_by_name(&graph, "B").expect("leaf B not found");
    let c = find_node_key_by_name(&graph, "C").expect("leaf C not found");

    let x = Array1::linspace(2000.0, 2010.0, 11);
    // means 2002, 2008, 2005 with precisions 1, 1, 2:
    // weighted mean = (2002 + 2008 + 2*2005) / 4 = 2005, which is a grid point.
    set_leaf_function(&graph, a, &x, gaussian_neglog(&x, 2002.0, 1.0))?;
    set_leaf_function(&graph, b, &x, gaussian_neglog(&x, 2008.0, 1.0))?;
    set_leaf_function(&graph, c, &x, gaussian_neglog(&x, 2005.0, 2.0))?;
    set_edge_branch_dist(&graph, a, 0.0);
    set_edge_branch_dist(&graph, b, 0.0);
    set_edge_branch_dist(&graph, c, 0.0);

    propagate_distributions_backward(&graph, None)?;

    let internal = find_node_key_by_name(&graph, "I").expect("internal node I not found");
    let node = graph.get_node(internal).expect("internal I exists");
    let payload = node.read_arc().payload().read_arc();
    let dist = payload
      .time_distribution()
      .as_ref()
      .expect("internal node should have a time distribution");
    let likely_time = dist.likely_time().expect("distribution should have a likely_time");

    // Mass re-windowing regrids the folded posterior, so 2005 need not be a stored grid point: the
    // discrete peak snaps to the nearest grid point, within one spacing of the analytic mean.
    let grid = dist.t();
    let spacing = grid[1] - grid[0];
    assert_abs_diff_eq!(likely_time, 2005.0, epsilon = spacing);

    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::graph::Graph;

    /// Give a node the date it was loaded with, and nothing else: the backward pass is what lifts
    /// it into the node's time distribution.
    pub(super) fn set_date_constraint(
      graph: &Graph<NodeTimetree, EdgeTimetree, ()>,
      key: GraphNodeKey,
      dist: Distribution<NegLog>,
    ) {
      let node = graph.get_node(key).expect("node exists");
      node.read_arc().payload().write_arc().date_constraint = Some(Arc::new(dist));
    }

    pub(super) fn set_leaf_time(graph: &Graph<NodeTimetree, EdgeTimetree, ()>, key: GraphNodeKey, time: f64) {
      let node = graph.get_node(key).expect("node exists");
      node.read_arc().payload().write_arc().time_distribution = Some(Arc::new(Distribution::point(time, 0.0)));
    }

    pub(super) fn set_edge_branch_dist(
      graph: &Graph<NodeTimetree, EdgeTimetree, ()>,
      target_key: GraphNodeKey,
      bl: f64,
    ) {
      for edge in graph.get_edges() {
        let edge_read = edge.read_arc();
        if edge_read.target() == target_key {
          edge_read
            .payload()
            .write_arc()
            .set_branch_length_distribution(Some(Arc::new(Distribution::point(bl, 0.0))));
        }
      }
    }

    /// A Gaussian likelihood stored in neg-log space is a parabola: `y(t) = precision*(t-mean)^2/2`.
    /// Its peak (the smallest ordinate) sits at `mean`, and a product of such messages is Gaussian
    /// with the precision-weighted mean.
    pub(super) fn gaussian_neglog(x: &Array1<f64>, mean: f64, precision: f64) -> Array1<f64> {
      x.mapv(|t| 0.5 * precision * (t - mean).powi(2))
    }

    /// Give a leaf a gridded `Function` time distribution so that its backward message is a
    /// `Function` (the operand kind that the common-grid fold must resample and add).
    pub(super) fn set_leaf_function(
      graph: &Graph<NodeTimetree, EdgeTimetree, ()>,
      key: GraphNodeKey,
      x: &Array1<f64>,
      y: Array1<f64>,
    ) -> Result<(), Report> {
      let dist = Distribution::function(x.clone(), y)?;
      let node = graph.get_node(key).expect("node exists");
      node.read_arc().payload().write_arc().time_distribution = Some(Arc::new(dist));
      Ok(())
    }

    pub(super) fn coalescent_model(tc: f64) -> Result<CoalescentModel, Report> {
      let lineage_counts = PiecewiseConstantFn::new(array![1900.0, 2100.0], array![1.0, 2.0, 0.0]);
      CoalescentModel::new(&lineage_counts, &Distribution::constant(tc))
    }
  }

  use helpers::*;
}
