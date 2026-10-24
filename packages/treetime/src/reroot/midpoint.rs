use crate::branch_lengths::branch_length_or_zero;
use crate::reroot::find_root::EdgePoint;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::iter::successors;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::split_fraction::SplitFraction;
use treetime_grid::grid::coordinates_coincide;
use treetime_utils::make_internal_report;

pub(crate) fn midpoint_root(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<Option<EdgePoint>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  match longest_path(&schedule, branch_lengths) {
    Some(path) if path.length > 0.0 => middle_of_path(&schedule, branch_lengths, &path).map(Some),
    _ => Ok(None),
  }
}

fn longest_path(
  schedule: &GraphSchedule<'_>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Option<LongestPath> {
  let mut deepest: BTreeMap<GraphNodeKey, Arm> = BTreeMap::new();
  let mut longest: Option<LongestPath> = None;
  for &key in schedule.postorder() {
    let arms = schedule
      .children(key)
      .iter()
      .map(|&(child_key, edge_key)| {
        let below = deepest[&child_key];
        Arm {
          length: below.length + branch_length_or_zero(branch_lengths, edge_key),
          leaf: below.leaf,
        }
      })
      .collect_vec();
    let is_root = schedule.parents(key).is_empty();
    let candidate = match arms.as_slice() {
      [] => None,
      [arm] => is_root.then_some(LongestPath {
        length: arm.length,
        end: arm.leaf,
      }),
      _ => two_longest_arms(&arms).map(|[first, second]| LongestPath {
        length: first.length + second.length,
        end: first.leaf,
      }),
    };
    if let Some(candidate) = candidate
      && longest.as_ref().is_none_or(|longest| candidate.length > longest.length)
    {
      longest = Some(candidate);
    }
    let arm = first_longest(arms).unwrap_or(Arm { length: 0.0, leaf: key });
    deepest.insert(key, arm);
  }
  longest
}

fn two_longest_arms(arms: &[Arm]) -> Option<[Arm; 2]> {
  let (first_index, first) = arms.iter().copied().enumerate().reduce(|best, current| {
    if current.1.length > best.1.length {
      current
    } else {
      best
    }
  })?;
  let others = arms
    .iter()
    .enumerate()
    .filter(|&(index, _)| index != first_index)
    .map(|(_, arm)| *arm);
  Some([first, first_longest(others)?])
}

fn first_longest(arms: impl IntoIterator<Item = Arm>) -> Option<Arm> {
  arms
    .into_iter()
    .reduce(|best, arm| if arm.length > best.length { arm } else { best })
}

fn middle_of_path(
  schedule: &GraphSchedule<'_>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  path: &LongestPath,
) -> Result<EdgePoint, Report> {
  let half = path.length / 2.0;
  let edges_up = successors(schedule.parents(path.end).first(), |&&(parent_key, _)| {
    schedule.parents(parent_key).first()
  });
  let mut covered = 0.0;
  for &(_, edge_key) in edges_up {
    let length = branch_length_or_zero(branch_lengths, edge_key);
    let reach = covered + length;
    let at_parent = coordinates_coincide(reach, half, path.length);
    if length > 0.0 && (at_parent || reach >= half) {
      let split = if at_parent {
        SplitFraction::PARENT_END
      } else {
        SplitFraction::new((reach - half) / length)?
      };
      return Ok(EdgePoint { edge: edge_key, split });
    }
    covered = reach;
  }
  Err(make_internal_report!(
    "The middle of the longest tip-to-tip path lies above the root"
  ))
}

#[derive(Clone, Copy, Debug)]
struct Arm {
  length: f64,
  leaf: GraphNodeKey,
}

struct LongestPath {
  length: f64,
  end: GraphNodeKey,
}
