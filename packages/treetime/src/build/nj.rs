use crate::build::query::Query;
use crate::gtr::jc_distance::jukes_cantor_distance;
use eyre::Report;
use itertools::Itertools;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};
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

/// Neighbor joining on a full distance matrix. With `top_hits = Some(m)`, each node compares itself only with
/// the `m` nodes that are best by the joining criterion when its list is made (FastTree's top hits), which replaces
/// the cubic scan of all pairs. Each node caches its best partner from its list. A joined node inherits the lists of
/// its children and enters the lists of its own partners, and all lists are rebuilt every `m` joins, so the search
/// costs $O(n^3 / m)$ instead of $O(n^3)$.
pub(crate) fn neighbor_joining(
  distances: &[f64],
  n_tips: usize,
  top_hits: Option<usize>,
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
  let mut alive = vec![true; n_tips];
  let mut sums = dist
    .chunks(n_tips.max(1))
    .map(|row| row.iter().sum::<f64>())
    .collect_vec();
  let mut hits = top_hits.map(|m| {
    (0..n_tips)
      .map(|row| best_partners(row, &active, active.len(), m, &dist, &sums))
      .collect_vec()
  });
  let mut best: Vec<Option<usize>> = vec![None; n_tips];
  while active.len() > 3 {
    let size = active.len().to_f64_exact();
    let q = |a: usize, b: usize| (size - 2.0) * dist[a * n_tips + b] - sums[a] - sums[b];
    let (_, first, second) = match (&mut hits, top_hits) {
      (Some(hits), Some(m)) => {
        let refresh = (n_tips - active.len()).is_multiple_of(m.max(1));
        let stale = active
          .iter()
          .copied()
          .filter(|&a| refresh || !hits[a].iter().any(|&b| alive[b]))
          .collect_vec();
        let fresh = stale
          .par_iter()
          .map(|&a| best_partners(a, &active, active.len(), m, &dist, &sums))
          .collect::<Vec<_>>();
        for (a, list) in stale.into_iter().zip(fresh) {
          hits[a] = list;
          best[a] = None;
        }
        for &a in &active {
          if best[a].is_none_or(|b| !alive[b]) {
            best[a] = hits[a]
              .iter()
              .copied()
              .filter(|&b| alive[b] && b != a)
              .min_by(|&x, &y| q(a, x).total_cmp(&q(a, y)).then(x.cmp(&y)));
          }
        }
        active
          .iter()
          .filter_map(|&a| best[a].map(|b| (q(a, b), a.min(b), a.max(b))))
          .min_by(|x, y| x.0.total_cmp(&y.0).then((x.1, x.2).cmp(&(y.1, y.2))))
          .ok_or_else(|| crate::make_internal_report!("Neighbor joining found no pair to join"))?
      },
      _ => active
        .iter()
        .tuple_combinations()
        .map(|(&a, &b)| (q(a, b), a, b))
        .min_by(|x, y| x.0.total_cmp(&y.0))
        .ok_or_else(|| crate::make_internal_report!("Neighbor joining found no pair to join"))?,
    };
    let joined = graph.add_node();
    graph.add_edge(joined, nodes[first])?;
    graph.add_edge(joined, nodes[second])?;
    alive[second] = false;
    active.retain(|&other| other != second);
    let mut joined_sum = 0.0;
    for &other in &active {
      if other != first {
        let value =
          0.5 * (dist[first * n_tips + other] + dist[second * n_tips + other] - dist[first * n_tips + second]);
        sums[other] += value - dist[first * n_tips + other] - dist[second * n_tips + other];
        joined_sum += value;
        dist[first * n_tips + other] = value;
        dist[other * n_tips + first] = value;
      }
    }
    sums[first] = joined_sum;
    if let (Some(hits), Some(m)) = (&mut hits, top_hits) {
      let candidates = hits[first]
        .iter()
        .chain(&hits[second])
        .copied()
        .filter(|&b| alive[b] && b != first)
        .sorted_unstable()
        .dedup()
        .collect_vec();
      hits[first] = best_partners(first, &candidates, active.len(), m, &dist, &sums);
      for &partner in &hits[first].clone() {
        if !hits[partner].contains(&first) {
          hits[partner].push(first);
        }
        best[partner] = None;
      }
      for &other in &active {
        if best[other].is_some_and(|b| b == first || b == second) {
          best[other] = None;
        }
      }
      hits[second] = vec![];
    }
    nodes[first] = joined;
  }
  let root = graph.add_node();
  for &remaining in &active {
    graph.add_edge(root, nodes[remaining])?;
  }
  graph.build()?;
  Ok((graph, samples))
}

fn best_partners(
  own: usize,
  candidates: &[usize],
  n_active: usize,
  m: usize,
  dist: &[f64],
  sums: &[f64],
) -> Vec<usize> {
  let n = sums.len();
  let scale = n_active.saturating_sub(2).to_f64_exact();
  let key = |other: usize| scale * dist[own * n + other] - sums[other];
  candidates
    .iter()
    .copied()
    .filter(|&other| other != own)
    .k_smallest_by(m, |&x, &y| key(x).total_cmp(&key(y)))
    .collect()
}
