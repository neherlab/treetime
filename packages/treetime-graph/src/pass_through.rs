use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SegmentTree<'a> {
  pub(crate) root: &'a str,
  pub(crate) children: BTreeMap<&'a str, Vec<&'a str>>,
  pub(crate) parent: BTreeMap<&'a str, (&'a str, Option<f64>)>,
}

impl<'a> SegmentTree<'a> {
  pub(crate) fn contains(&self, name: &str) -> bool {
    self.root == name || self.parent.contains_key(name)
  }

  pub(crate) fn preorder(&self) -> Vec<&'a str> {
    let mut order = Vec::with_capacity(self.parent.len() + 1);
    let mut stack = vec![self.root];
    while let Some(name) = stack.pop() {
      order.push(name);
      if let Some(children) = self.children.get(name) {
        stack.extend(children.iter().rev());
      }
    }
    order
  }

  fn path_between(&self, ancestor: &str, descendant: &str) -> Option<Vec<(&'a str, Option<f64>)>> {
    let mut steps = vec![];
    let mut node = descendant;
    loop {
      let &(parent, length) = self.parent.get(node)?;
      steps.push((parent, length));
      if parent == ancestor {
        steps.reverse();
        return Some(steps);
      }
      node = parent;
    }
  }

  fn reroute(&mut self, ancestor: &'a str, descendant: &'a str, path: &[&'a str], lengths: &[Option<f64>]) {
    let chain = std::iter::once(ancestor)
      .chain(path.iter().copied())
      .chain([descendant])
      .collect_vec();
    if let Some(children) = self.children.get_mut(ancestor)
      && let Some(slot) = children.iter_mut().find(|child| **child == descendant)
    {
      *slot = chain[1];
    }
    for (index, (&parent, &child)) in chain.iter().tuple_windows().enumerate() {
      self
        .parent
        .insert(child, (parent, lengths.get(index).copied().flatten()));
      if parent != ancestor {
        self.children.insert(parent, vec![child]);
      }
    }
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PassThroughReport {
  pub inserted_nodes: usize,
  pub ambiguous_edges: usize,
  pub conflicting_edges: usize,
}

pub(crate) fn insert_pass_through_nodes(trees: &mut [SegmentTree<'_>]) -> PassThroughReport {
  let mut report = PassThroughReport::default();
  let mut ambiguous = BTreeSet::new();
  let mut conflicting = BTreeSet::new();
  loop {
    let mut changed = false;
    for segment in 0..trees.len() {
      let edges = trees[segment]
        .parent
        .iter()
        .map(|(&child, &(parent, length))| (parent, child, length))
        .collect_vec();
      for (parent, child, length) in edges {
        let candidates = trees
          .iter()
          .enumerate()
          .filter(|&(other, _)| other != segment)
          .filter_map(|(_, tree)| tree.path_between(parent, child))
          .filter(|path| path.len() > 1)
          .collect_vec();
        let Some(longest) = candidates.iter().max_by_key(|path| path.len()) else {
          continue;
        };
        let longest_nodes = interior(longest);
        if !candidates
          .iter()
          .all(|path| is_subsequence(&interior(path), &longest_nodes))
        {
          ambiguous.insert((segment, parent, child));
          continue;
        }
        if longest_nodes.iter().any(|name| trees[segment].contains(name)) {
          conflicting.insert((segment, parent, child));
          continue;
        }
        let steps = longest.iter().map(|&(_, length)| length).collect_vec();
        let lengths = split_length(length, &steps);
        trees[segment].reroute(parent, child, &longest_nodes, &lengths);
        report.inserted_nodes += longest_nodes.len();
        changed = true;
      }
    }
    if !changed {
      break;
    }
  }
  report.ambiguous_edges = ambiguous.len();
  report.conflicting_edges = conflicting.len();
  report
}

fn interior<'a>(path: &[(&'a str, Option<f64>)]) -> Vec<&'a str> {
  path.iter().skip(1).map(|&(name, _)| name).collect()
}

fn is_subsequence(short: &[&str], long: &[&str]) -> bool {
  let mut long = long.iter();
  short.iter().all(|name| long.any(|other| other == name))
}

fn split_length(total: Option<f64>, steps: &[Option<f64>]) -> Vec<Option<f64>> {
  let Some(total) = total else {
    return vec![None; steps.len()];
  };
  let known = steps.iter().copied().collect::<Option<Vec<f64>>>();
  let sum = known.as_ref().map_or(0.0, |lengths| lengths.iter().sum::<f64>());
  match known {
    Some(lengths) if sum > 0.0 => lengths.iter().map(|length| Some(total * length / sum)).collect(),
    _ => {
      let count = f64::from(u32::try_from(steps.len()).unwrap_or(u32::MAX));
      vec![Some(total / count); steps.len()]
    },
  }
}
