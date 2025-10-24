use crate::distribution::distribution::Distribution;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::GraphNodeBackward;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::Report;
use ndarray::Array1;
use std::sync::Arc;

const PARENT_GRID_SIZE: usize = 200;
const RANGE_SAMPLE_POINTS: usize = 64;
const EPS: f64 = 1e-9;

/// Propagates time distributions backward from leaves to root.
pub fn propagate_distributions_backward(graph: &GraphAncestral) -> Result<(), Report> {
  graph.par_iter_breadth_first_backward(|mut node| {
    propagate_distributions_backward_single_node(&mut node).unwrap();
    GraphTraversalContinuation::Continue
  });
  Ok(())
}

/// Computes time distribution for a single internal node from its children.
fn propagate_distributions_backward_single_node(
  node: &mut GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  if node.is_leaf {
    return Ok(());
  }

  let mut result: Option<Distribution> = None;
  for (child, edge) in &node.children {
    let child = child.read_arc();
    let mut edge = edge.write_arc();

    if let (Some(branch_dist), Some(child_time_dist)) = (&edge.branch_length_distribution, &child.time_distribution) {
      // Compute parent time distribution from child time and branch length
      let new = compute_parent_message(child_time_dist.as_ref(), branch_dist.as_ref())?;
      let new_arc = Arc::new(new);

      // Store message on edge
      edge.msg_to_parent = Some(Arc::clone(&new_arc));

      // Multiply messages from all children to find intersection of constraints
      result = Some(if let Some(current) = &result {
        multiply_parent_distributions(current, &new_arc)?
      } else {
        (*new_arc).clone()
      });
    }
  }

  // Store final distribution on node
  if let Some(dist) = result {
    node.payload.time_distribution = Some(Arc::new(dist));
  }

  Ok(())
}

/// Computes parent time distribution given child time and branch length: t_parent = t_child - d.
fn compute_parent_message(child: &Distribution, branch: &Distribution) -> Result<Distribution, Report> {
  if matches!(child, Distribution::Empty) || matches!(branch, Distribution::Empty) {
    return Ok(Distribution::empty());
  }

  // Special case: child is point mass, directly shift branch distribution
  if let Distribution::Point(point) = child {
    return shift_branch_distribution(point.t(), point.amplitude(), branch);
  }

  // Get support bounds for child distribution
  let (child_min, child_max) = match distribution_support_bounds(child) {
    (Some(min), Some(max)) => (min, max),
    _ => return Ok(Distribution::empty()),
  };

  // Sample branch distribution into discrete points
  let branch_points = sample_distribution_points(branch, RANGE_SAMPLE_POINTS);
  if branch_points.is_empty() {
    return Ok(Distribution::empty());
  }

  // Find branch duration range
  let branch_min = branch_points.iter().map(|(d, _)| *d).fold(f64::INFINITY, f64::min);
  let branch_max = branch_points.iter().map(|(d, _)| *d).fold(f64::NEG_INFINITY, f64::max);

  // Compute parent time range: t_parent = t_child - duration
  let parent_min = child_min - branch_max;
  let parent_max = child_max - branch_min;

  // Special case: parent range collapses to point
  if (parent_max - parent_min).abs() <= EPS {
    let amplitude = branch_points
      .iter()
      .map(|(_, w)| *w)
      .fold(0.0, |acc, w| acc + w)
      .max(0.0);
    return Ok(Distribution::point(parent_min, amplitude));
  }

  // Create uniform grid over parent time range
  let step = (parent_max - parent_min) / (PARENT_GRID_SIZE - 1) as f64;
  let parent_times: Vec<f64> = (0..PARENT_GRID_SIZE).map(|i| parent_min + step * i as f64).collect();

  // Accumulate probabilities in log-space for numerical stability
  let mut parent_log = vec![f64::NEG_INFINITY; parent_times.len()];

  // Marginalize over branch durations: P(t_parent) = sum_d P(d) * P(t_child = t_parent + d)
  for (duration, weight) in branch_points.iter() {
    if *weight <= 0.0 {
      continue;
    }
    let log_w = weight.ln();

    for (i, &t_parent) in parent_times.iter().enumerate() {
      // Compute implied child time
      let child_time = t_parent + *duration;
      let value = distribution_value(child, child_time);
      if value <= 0.0 {
        continue;
      }

      // Accumulate in log-space
      let total_log = log_w + value.ln();
      parent_log[i] = log_add(parent_log[i], total_log);
    }
  }

  if parent_log.iter().all(|v| !v.is_finite()) {
    return Ok(Distribution::empty());
  }

  // Normalize and convert back to linear space
  let max_log = parent_log
    .iter()
    .copied()
    .filter(|v| v.is_finite())
    .fold(f64::NEG_INFINITY, f64::max);

  let values: Vec<f64> = parent_log
    .iter()
    .map(|&log_v| {
      if log_v.is_finite() {
        (log_v - max_log).exp()
      } else {
        0.0
      }
    })
    .collect();

  Distribution::function(Array1::from_vec(parent_times), Array1::from_vec(values))
}

/// Shifts branch distribution backward when child is point mass: t_parent = t_child - d.
fn shift_branch_distribution(
  child_time: f64,
  child_weight: f64,
  branch: &Distribution,
) -> Result<Distribution, Report> {
  let points = sample_distribution_points(branch, RANGE_SAMPLE_POINTS);
  if points.is_empty() {
    return Ok(Distribution::empty());
  }

  // Transform each (duration, probability) to (parent_time, probability)
  let mut transformed: Vec<(f64, f64)> = points
    .into_iter()
    .map(|(d, w)| (child_time - d, w * child_weight))
    .collect();

  // Sort by time to create valid distribution
  transformed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

  let times: Vec<f64> = transformed.iter().map(|(t, _)| *t).collect();
  let values: Vec<f64> = transformed.iter().map(|(_, v)| *v).collect();
  Distribution::function(Array1::from_vec(times), Array1::from_vec(values))
}

/// Multiplies two distributions to find intersection of time constraints.
fn multiply_parent_distributions(a: &Distribution, b: &Distribution) -> Result<Distribution, Report> {
  if matches!(a, Distribution::Empty) || matches!(b, Distribution::Empty) {
    return Ok(Distribution::empty());
  }

  // Find support bounds for both distributions
  let (min_a, max_a) = distribution_support_bounds(a);
  let (min_b, max_b) = distribution_support_bounds(b);
  let (Some(min_a), Some(max_a)) = (min_a, max_a) else {
    return Ok(Distribution::empty());
  };
  let (Some(min_b), Some(max_b)) = (min_b, max_b) else {
    return Ok(Distribution::empty());
  };

  // Compute overlapping region
  let overlap_min = min_a.max(min_b);
  let overlap_max = max_a.min(max_b);
  if overlap_min >= overlap_max {
    return Ok(Distribution::empty());
  }

  // Special case: overlap is single point
  if (overlap_max - overlap_min).abs() <= EPS {
    let value = distribution_value(a, overlap_min) * distribution_value(b, overlap_min);
    if value <= 0.0 {
      return Ok(Distribution::empty());
    }
    return Ok(Distribution::point(overlap_min, value));
  }

  // Create grid over overlap region
  let step = (overlap_max - overlap_min) / (PARENT_GRID_SIZE - 1) as f64;
  let times: Vec<f64> = (0..PARENT_GRID_SIZE).map(|i| overlap_min + step * i as f64).collect();
  let mut log_values = Vec::with_capacity(times.len());

  // Multiply distributions pointwise in log-space
  for &t in &times {
    let va = distribution_value(a, t);
    let vb = distribution_value(b, t);
    if va <= 0.0 || vb <= 0.0 {
      log_values.push(f64::NEG_INFINITY);
    } else {
      log_values.push(va.ln() + vb.ln());
    }
  }

  if log_values.iter().all(|v| !v.is_finite()) {
    return Ok(Distribution::empty());
  }

  // Normalize and convert to linear space
  let max_log = log_values
    .iter()
    .copied()
    .filter(|v| v.is_finite())
    .fold(f64::NEG_INFINITY, f64::max);

  let values: Vec<f64> = log_values
    .iter()
    .map(|&log_v| {
      if log_v.is_finite() {
        (log_v - max_log).exp()
      } else {
        0.0
      }
    })
    .collect();

  Distribution::function(Array1::from_vec(times), Array1::from_vec(values))
}

/// Samples distribution into discrete (time, probability) points.
fn sample_distribution_points(dist: &Distribution, samples: usize) -> Vec<(f64, f64)> {
  match dist {
    Distribution::Empty => vec![],
    Distribution::Point(point) => vec![(point.t(), point.amplitude())],
    Distribution::Range(range) => {
      if (range.end() - range.start()).abs() <= EPS {
        vec![(range.start(), range.amplitude())]
      } else {
        let n = samples.max(2);
        (0..n)
          .map(|i| {
            let ratio = i as f64 / (n - 1) as f64;
            (range.start() + (range.end() - range.start()) * ratio, range.amplitude())
          })
          .collect()
      }
    },
    Distribution::Function(function) => function
      .t()
      .iter()
      .zip(function.y().iter())
      .map(|(&t, &y)| (t, y))
      .collect(),
  }
}

/// Evaluates distribution probability at given time.
fn distribution_value(dist: &Distribution, time: f64) -> f64 {
  match dist {
    Distribution::Empty => 0.0,
    Distribution::Point(point) => {
      if (time - point.t()).abs() <= 1e-6 {
        point.amplitude()
      } else {
        0.0
      }
    },
    Distribution::Range(range) => {
      if time >= range.start() - EPS && time <= range.end() + EPS {
        range.amplitude()
      } else {
        0.0
      }
    },
    Distribution::Function(function) => function.interp(time).unwrap_or(0.0),
  }
}

/// Returns (min, max) time bounds where distribution has support.
fn distribution_support_bounds(dist: &Distribution) -> (Option<f64>, Option<f64>) {
  match dist {
    Distribution::Empty => (None, None),
    Distribution::Point(point) => {
      let t = point.t();
      (Some(t), Some(t))
    },
    Distribution::Range(range) => (Some(range.start()), Some(range.end())),
    Distribution::Function(function) => {
      if function.t().is_empty() {
        (None, None)
      } else {
        let t = function.t();
        (Some(t[0]), Some(t[t.len() - 1]))
      }
    },
  }
}

/// Computes log(exp(a) + exp(b)) numerically stable (log-sum-exp trick).
fn log_add(a: f64, b: f64) -> f64 {
  if !a.is_finite() {
    return b;
  }
  if !b.is_finite() {
    return a;
  }
  let max = a.max(b);
  max + ((a - max).exp() + (b - max).exp()).ln()
}

/// Estimates average spacing between consecutive time points.
fn estimate_step(points: &[(f64, f64)]) -> Option<f64> {
  if points.len() < 2 {
    return None;
  }

  let mut deltas = Vec::new();
  for window in points.windows(2) {
    let delta = (window[1].0 - window[0].0).abs();
    if delta > 0.0 {
      deltas.push(delta);
    }
  }

  if deltas.is_empty() {
    None
  } else {
    Some(deltas.iter().copied().sum::<f64>() / deltas.len() as f64)
  }
}
