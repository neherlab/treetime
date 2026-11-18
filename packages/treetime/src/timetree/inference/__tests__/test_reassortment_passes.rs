#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::optimize::indel::IndelTerm;
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::pretty_assert_ulps_eq;
  use crate::progress::NoopProgress;
  use crate::test_utils::{
    RecordingLog, find_node_key_by_name, full_posterior, point_date_constraints, unknown_branches,
  };
  use crate::timetree::inference::backward_pass::{BackwardInputs, ParentMessages};
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::branch_length_likelihood::compute_branch_length_distribution;
  use crate::timetree::inference::forward_pass::{ForwardPassInputs, propagate_distributions_forward};
  use crate::timetree::inference::grid_request::GridOrigin;
  use crate::timetree::inference::result::{BranchLikelihood, NodeBackward};
  use crate::timetree::inference::runner::{GRID_POINTS, sweep_time_passes};
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::Array2;
  use pretty_assertions::assert_eq;
  use smallvec::smallvec;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::Arc;
  use treetime_distribution::{Distribution, MassProfile, NegLog, TailMass};
  use treetime_graph::key_map::{EdgeMap, NodeMap};
  use treetime_graph::merge_by_name::{MergedGraph, NamedTree, merge_trees_by_name};
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::TREE_DIALECT_DEFAULT;
  use treetime_io::tree::tree_read;
  use treetime_ops::FftPlans;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::make_report;

  const RATE: f64 = 2.0;
  const SITES: usize = 1000;
  const TIPS: [(&str, f64); 5] = [
    ("A", 2012.0),
    ("B", 2011.0),
    ("D", 2012.5),
    ("E", 2012.0),
    ("F", 2009.0),
  ];
  const POLYTREE: [&str; 2] = ["(A:1,(D:1,E:1)C:1)P1;", "(B:1,(D:1,E:1)C:1)P2;"];
  const LOOP: [&str; 2] = ["((A:1,(D:1,E:1)C:1)P1:1,F:1)R;", "((B:1,(D:1,E:1)C:1)P2:1,F:1)R;"];
  const COUNTS: [((&str, &str), usize); 9] = [
    (("P1", "A"), 4),
    (("P1", "C"), 6),
    (("P2", "B"), 3),
    (("P2", "C"), 1),
    (("C", "D"), 2),
    (("C", "E"), 3),
    (("R", "P1"), 2),
    (("R", "P2"), 2),
    (("R", "F"), 8),
  ];

  #[test]
  fn test_reassortment_passes_on_a_graph_without_cycles_match_the_exact_posterior() -> Result<(), Report> {
    let grid = grid(2000.0, 0.02, 651);
    let exact = exact_polytree(&grid);
    let all = passes_cdfs(&POLYTREE, ParentMessages::AllEdges, &grid, &["P1", "P2", "C"])?;
    let children = passes_cdfs(&POLYTREE, ParentMessages::ChildrenOnly, &grid, &["P1", "P2", "C"])?;
    let errors =
      |cdfs: &BTreeMap<String, Vec<f64>>| ["P1", "P2", "C"].map(|node| max_difference(&cdf(&exact[node]), &cdfs[node]));
    let (all, children) = (errors(&all), errors(&children));
    assert!(all.iter().all(|&error| error < 5e-3), "all edges: {all:?}");
    assert!(
      children[0] > 4.0 * all[0],
      "P1 without the other parent: {children:?} against {all:?}"
    );
    Ok(())
  }

  #[test]
  fn test_reassortment_passes_on_a_graph_with_a_cycle_stay_close_to_the_exact_posterior_with_the_messages_of_both_parents()
  -> Result<(), Report> {
    let grid = grid(1995.0, 0.075, 241);
    let exact = exact_loop(&grid);
    let nodes = ["R", "P1", "P2", "C"];
    let all = passes_cdfs(&LOOP, ParentMessages::AllEdges, &grid, &nodes)?;
    let children = passes_cdfs(&LOOP, ParentMessages::ChildrenOnly, &grid, &nodes)?;
    let errors = |cdfs: &BTreeMap<String, Vec<f64>>| nodes.map(|node| max_difference(&cdf(&exact[node]), &cdfs[node]));
    let (all, children) = (errors(&all), errors(&children));
    assert!(all.iter().all(|&error| error < 0.015), "all edges: {all:?}");
    assert!(
      children[1] > 3.0 * all[1] && children[2] > 3.0 * all[2],
      "the parents without each other's message: {children:?} against {all:?}"
    );
    Ok(())
  }

  #[test]
  fn test_reassortment_passes_date_a_node_from_its_subtree_when_its_parents_leave_no_common_time() -> Result<(), Report>
  {
    let MergedGraph { graph, names, .. } = merge(&POLYTREE)?;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let mut backward: NodeMap<NodeBackward> = graph
      .get_nodes()
      .map(|node| (node.key(), NodeBackward::default()))
      .collect();
    for (name, time) in [("P1", 2000.0), ("P2", 2010.0), ("C", 2015.0)] {
      backward.entry(key(name)).or_default().subtree = Some(Arc::new(Distribution::point(time, 0.0)));
    }
    let child = backward.entry(key("C")).or_default();
    child.to_parents = smallvec![None, None];
    child.to_parent_origins = smallvec![GridOrigin::default(), GridOrigin::default()];
    let mut branches = unknown_branches(&graph);
    for edge in graph.get_edges().filter(|edge| &names[edge.target()] == "C") {
      let branch = BranchLikelihood {
        distribution: Some(Arc::new(Distribution::point(1.0, 0.0))),
        time_length: Some(1.0),
      };
      branches.insert(edge.key(), branch);
    }
    let constraints = DateConstraints::default();
    let bad_branches = graph.get_nodes().map(|node| (node.key(), false)).collect();
    let fft_plans = FftPlans::new();
    let inputs = ForwardPassInputs {
      constraints: &constraints,
      coalescent: None,
      bad_branches: &bad_branches,
      names: &names,
      branches: &branches,
      max_points: MaxGridPoints::default(),
      fft_plans: &fft_plans,
    };
    let log = RecordingLog::default();

    let posterior = propagate_distributions_forward(&graph, &inputs, backward, full_posterior, &log)?;

    let expected_warnings = vec![
      "Timetree forward pass: 1 node(s) with several parents received messages from their parents and their \
       subtree that leave no time with any probability, so each took its time from its subtree alone, or stayed \
       undated when nothing below it is dated. The parents place such a node at incompatible times, for example a \
       reassortant whose segment trees disagree on its age: C"
        .to_owned(),
    ];
    assert_eq!(expected_warnings, log.warnings());
    let time = posterior[&key("C")]
      .time
      .ok_or_else(|| make_report!("node C has no time"))?;
    pretty_assert_ulps_eq!(2015.0, time, max_ulps = 4);
    Ok(())
  }

  fn grid(start: f64, step: f64, len: usize) -> Vec<f64> {
    (0..len).map(|i| start + step * i.to_f64_exact()).collect()
  }

  fn count(parent: &str, child: &str) -> usize {
    COUNTS
      .iter()
      .find(|((p, c), _)| *p == parent && *c == child)
      .map_or(0, |(_, k)| *k)
  }

  fn tip(name: &str) -> f64 {
    TIPS
      .iter()
      .find(|(tip, _)| *tip == name)
      .map_or(f64::NAN, |(_, date)| *date)
  }

  fn branch_lh(gtr: &GTR, tau: f64, k: usize) -> f64 {
    if tau < 0.0 || (tau == 0.0 && k > 0) {
      return 0.0;
    }
    let b = RATE / SITES.to_f64_exact() * tau;
    let eigvals = gtr.rate_eigvals();
    let p = |a: usize, c: usize| {
      (0..eigvals.len())
        .map(|j| gtr.v()[[a, j]] * gtr.v_inv()[[j, c]] * (eigvals[j] * b).exp())
        .sum::<f64>()
    };
    let (same, other) = (p(0, 0), p(0, 1));
    let differing = if k == 0 { 0.0 } else { k.to_f64_exact() * other.ln() };
    ((SITES - k).to_f64_exact() * same.ln() + differing).exp()
  }

  fn factor<'g>(gtr: &'g GTR, parent: &str, child: &str) -> impl Fn(f64, f64) -> f64 + 'g {
    let k = count(parent, child);
    move |t_parent, t_child| branch_lh(gtr, t_child - t_parent, k)
  }

  fn exact_polytree(grid: &[f64]) -> BTreeMap<&'static str, Vec<f64>> {
    let gtr = jc69(JC69Params::default()).expect("JC69 builds");
    let (p1c, p2c) = (factor(&gtr, "P1", "C"), factor(&gtr, "P2", "C"));
    let tip_a = |t: f64| factor(&gtr, "P1", "A")(t, tip("A"));
    let tip_b = |t: f64| factor(&gtr, "P2", "B")(t, tip("B"));
    let tips_c = |t: f64| factor(&gtr, "C", "D")(t, tip("D")) * factor(&gtr, "C", "E")(t, tip("E"));
    let from_p1: Vec<f64> = grid
      .iter()
      .map(|&c| grid.iter().map(|&p| tip_a(p) * p1c(p, c)).sum())
      .collect();
    let from_p2: Vec<f64> = grid
      .iter()
      .map(|&c| grid.iter().map(|&p| tip_b(p) * p2c(p, c)).sum())
      .collect();
    let marginal_c: Vec<f64> = grid
      .iter()
      .enumerate()
      .map(|(i, &t)| tips_c(t) * from_p1[i] * from_p2[i])
      .collect();
    let p1 = grid
      .iter()
      .map(|&p| {
        tip_a(p)
          * grid
            .iter()
            .enumerate()
            .map(|(i, &t)| p1c(p, t) * tips_c(t) * from_p2[i])
            .sum::<f64>()
      })
      .collect();
    let p2 = grid
      .iter()
      .map(|&p| {
        tip_b(p)
          * grid
            .iter()
            .enumerate()
            .map(|(i, &t)| p2c(p, t) * tips_c(t) * from_p1[i])
            .sum::<f64>()
      })
      .collect();
    BTreeMap::from([("P1", p1), ("P2", p2), ("C", marginal_c)])
  }

  fn exact_loop(grid: &[f64]) -> BTreeMap<&'static str, Vec<f64>> {
    let gtr = jc69(JC69Params::default()).expect("JC69 builds");
    let size = grid.len();
    let at = |f: &dyn Fn(f64, f64) -> f64| Array2::from_shape_fn((size, size), |(i, j)| f(grid[i], grid[j]));
    let (rp1, rp2) = (at(&factor(&gtr, "R", "P1")), at(&factor(&gtr, "R", "P2")));
    let (p1c, p2c) = (at(&factor(&gtr, "P1", "C")), at(&factor(&gtr, "P2", "C")));
    let r_tip: Vec<f64> = grid.iter().map(|&t| factor(&gtr, "R", "F")(t, tip("F"))).collect();
    let tip_a: Vec<f64> = grid.iter().map(|&t| factor(&gtr, "P1", "A")(t, tip("A"))).collect();
    let tip_b: Vec<f64> = grid.iter().map(|&t| factor(&gtr, "P2", "B")(t, tip("B"))).collect();
    let tips_c: Vec<f64> = grid
      .iter()
      .map(|&t| factor(&gtr, "C", "D")(t, tip("D")) * factor(&gtr, "C", "E")(t, tip("E")))
      .collect();
    let (mut marginal_r, mut p1, mut p2, mut marginal_c) =
      (vec![0.0; size], vec![0.0; size], vec![0.0; size], vec![0.0; size]);
    for i in 0..size {
      for j in 0..size {
        let below: f64 = (0..size).map(|k| p1c[[i, k]] * p2c[[j, k]] * tips_c[k]).sum();
        let above: f64 = (0..size).map(|q| r_tip[q] * rp1[[q, i]] * rp2[[q, j]]).sum();
        let weight = tip_a[i] * tip_b[j];
        p1[i] += weight * below * above;
        p2[j] += weight * below * above;
        for q in 0..size {
          marginal_r[q] += weight * below * r_tip[q] * rp1[[q, i]] * rp2[[q, j]];
        }
        for k in 0..size {
          marginal_c[k] += weight * above * p1c[[i, k]] * p2c[[j, k]] * tips_c[k];
        }
      }
    }
    BTreeMap::from([("R", marginal_r), ("P1", p1), ("P2", p2), ("C", marginal_c)])
  }

  fn merge(trees: &[&str]) -> Result<MergedGraph, Report> {
    let parsed = trees
      .iter()
      .map(|tree| tree_read(tree.as_bytes(), TREE_DIALECT_DEFAULT))
      .collect::<Result<Vec<_>, _>>()?;
    let labels = (1..=parsed.len()).map(|index| format!("tree {index}")).collect_vec();
    let named = parsed
      .iter()
      .zip(&labels)
      .map(|(tree, label)| NamedTree {
        label,
        graph: &tree.graph,
        names: &tree.names,
        unnamed_nodes: &tree.unnamed_nodes,
        branch_lengths: &tree.branch_lengths,
      })
      .collect_vec();
    merge_trees_by_name(&named)
  }

  fn branch_distribution(gtr: &GTR, k: usize) -> Result<Arc<Distribution<NegLog>>, Report> {
    let n_states = gtr.v().nrows();
    let mut coefficients = Array2::<f64>::zeros((SITES, n_states));
    for site in 0..SITES {
      let child = usize::from(site < k);
      for j in 0..n_states {
        coefficients[[site, j]] = gtr.v()[[0, j]] * gtr.v_inv()[[j, child]];
      }
    }
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));
    let one_mutation = 1.0 / SITES.to_f64_exact();
    compute_branch_length_distribution(
      &[contribution],
      IndelTerm::new(0, 0.0),
      (k.to_f64_exact() * one_mutation).max(one_mutation),
      one_mutation,
      GRID_POINTS,
      RATE * one_mutation,
      MaxGridPoints::default(),
    )
  }

  fn passes_cdfs(
    trees: &[&str],
    parent_messages: ParentMessages,
    grid: &[f64],
    nodes: &[&str],
  ) -> Result<BTreeMap<String, Vec<f64>>, Report> {
    let gtr = jc69(JC69Params::default())?;
    let merged = merge(trees)?;
    let MergedGraph { graph, names, .. } = merged;
    let branches: EdgeMap<BranchLikelihood> = graph
      .get_edges()
      .map(|edge| {
        let k = count(&names[edge.source()], &names[edge.target()]);
        let distribution = branch_distribution(&gtr, k)?;
        let time_length = distribution.likely_time()?;
        Ok((
          edge.key(),
          BranchLikelihood {
            distribution: Some(distribution),
            time_length,
          },
        ))
      })
      .collect::<Result<_, Report>>()?;
    let dated = TIPS
      .iter()
      .copied()
      .filter(|(name, _)| graph.get_nodes().any(|node| &names[node.key()] == *name))
      .collect_vec();
    let constraints = point_date_constraints(&graph, &names, &dated);
    let bad_branches = derive_bad_branches(
      &graph,
      &constraints,
      &bad_leaves(&graph, &constraints, &BTreeSet::new()),
    )?;
    let fft_plans = FftPlans::new();
    let backward = BackwardInputs {
      constraints: &constraints,
      coalescent_model: None,
      bad_branches: &bad_branches,
      branches: &branches,
      max_points: MaxGridPoints::default(),
      fft_plans: &fft_plans,
    };
    let forward = ForwardPassInputs {
      constraints: &constraints,
      coalescent: None,
      bad_branches: &bad_branches,
      names: &names,
      branches: &branches,
      max_points: MaxGridPoints::default(),
      fft_plans: &fft_plans,
    };
    let (_, posterior) = sweep_time_passes(
      &graph,
      &backward,
      &forward,
      parent_messages,
      full_posterior,
      |posterior| posterior.likely_time,
      &NoopProgress,
    )?;
    nodes
      .iter()
      .map(|&name| {
        let key = graph
          .get_nodes()
          .map(|node| node.key())
          .find(|&key| &names[key] == name)
          .ok_or_else(|| make_report!("node {name} missing"))?;
        let distribution = posterior[&key]
          .distribution
          .as_deref()
          .ok_or_else(|| make_report!("node {name} has no posterior"))?;
        Ok((name.to_owned(), posterior_cdf(distribution, grid)?))
      })
      .collect()
  }

  fn posterior_cdf(distribution: &Distribution<NegLog>, grid: &[f64]) -> Result<Vec<f64>, Report> {
    let Distribution::Function(f) = distribution else {
      return Err(make_report!("an internal node posterior must be a function"));
    };
    let profile = MassProfile::new(f)?;
    let fractions = profile
      .fractions(TailMass::Required)
      .ok_or_else(|| make_report!("an internal node posterior must have a finite mass"))?;
    Ok(grid.iter().map(|&t| fractions.cdf(t)).collect())
  }

  fn cdf(density: &[f64]) -> Vec<f64> {
    let mut cumulative = vec![0.0; density.len()];
    for i in 1..density.len() {
      cumulative[i] = cumulative[i - 1] + f64::midpoint(density[i], density[i - 1]);
    }
    let total = cumulative.last().copied().unwrap_or(1.0);
    cumulative.iter().map(|c| c / total).collect()
  }

  fn max_difference(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max)
  }
}
