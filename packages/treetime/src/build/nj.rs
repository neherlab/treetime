use crate::build::query::Query;
use crate::gtr::jc_distance::jukes_cantor_distance;
use eyre::Report;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_utils::cast::ToF64Exact;

pub(crate) fn jc_distances(queries: &[Query], reference: &Seq, n_states: usize) -> Vec<f64> {
  let n = queries.len();
  let length = reference.len();
  (0..n * n)
    .into_par_iter()
    .map(|index| {
      let (i, j) = (index.div_euclid(n), index.rem_euclid(n));
      if i == j {
        return 0.0;
      }
      let compared = length
        .saturating_sub(queries[i].n_missing() + queries[j].n_missing())
        .max(1);
      let p = queries[i].distance(&queries[j], reference).to_f64_exact() / compared.to_f64_exact();
      jukes_cantor_distance(p, n_states)
    })
    .collect()
}

pub(crate) fn neighbor_joining(
  distances: &[f64],
  n_tips: usize,
) -> Result<(Graph, BTreeMap<GraphNodeKey, usize>), Report> {
  let mut graph = Graph::new();
  let mut samples = BTreeMap::new();
  let mut nodes = (0..n_tips)
    .map(|sample| {
      let node = graph.add_node();
      samples.insert(node, sample);
      node
    })
    .collect_vec();
  let mut dist = distances.to_vec();
  let mut active = (0..n_tips).collect_vec();
  while active.len() > 3 {
    let size = active.len().to_f64_exact();
    let sums = active
      .iter()
      .map(|&row| active.iter().map(|&col| dist[row * n_tips + col]).sum::<f64>())
      .collect_vec();
    let mut best = (f64::INFINITY, 0, 1);
    for left in 0..active.len() {
      for right in left + 1..active.len() {
        let score = (size - 2.0) * dist[active[left] * n_tips + active[right]] - sums[left] - sums[right];
        if score < best.0 {
          best = (score, left, right);
        }
      }
    }
    let (_, left, right) = best;
    let (first, second) = (active[left], active[right]);
    let joined = graph.add_node();
    graph.add_edge(joined, nodes[first])?;
    graph.add_edge(joined, nodes[second])?;
    for &other in &active {
      if other != first && other != second {
        let value =
          0.5 * (dist[first * n_tips + other] + dist[second * n_tips + other] - dist[first * n_tips + second]);
        dist[first * n_tips + other] = value;
        dist[other * n_tips + first] = value;
      }
    }
    nodes[first] = joined;
    active.remove(right);
  }
  let root = graph.add_node();
  for &remaining in &active {
    graph.add_edge(root, nodes[remaining])?;
  }
  graph.build()?;
  Ok((graph, samples))
}
