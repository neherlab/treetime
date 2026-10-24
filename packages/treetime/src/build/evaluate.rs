use crate::alphabet::alphabet::Alphabet;
use crate::make_error;
use crate::partition::fitch::passes::create_fitch_partition;
use crate::seq::alignment::NodeSeqInput;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::clades::{TipIds, graph_clade_keys};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::Seq;

pub fn fitch_edge_counts(
  graph: &Graph,
  names: &NodeNames,
  seqs: &BTreeMap<&str, &Seq>,
  alphabet: &Alphabet,
) -> Result<BTreeMap<GraphEdgeKey, usize>, Report> {
  Ok(
    fitch_edge_subs(graph, names, seqs, alphabet)?
      .into_iter()
      .map(|(edge, subs)| (edge, subs.len()))
      .collect(),
  )
}

pub fn fitch_edge_subs(
  graph: &Graph,
  names: &NodeNames,
  seqs: &BTreeMap<&str, &Seq>,
  alphabet: &Alphabet,
) -> Result<BTreeMap<GraphEdgeKey, Vec<Sub>>, Report> {
  let inputs = graph
    .get_leaves()
    .map(|leaf| {
      let name = &names[leaf.key()];
      let Some(&seq) = seqs.get(name) else {
        return make_error!("Tree tip '{name}' has no sequence");
      };
      Ok((
        leaf.key(),
        NodeSeqInput {
          name: name.to_owned(),
          seq: Some(seq.clone()),
        },
      ))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let fitch = create_fitch_partition(graph, alphabet.clone(), inputs)?;
  Ok(
    fitch
      .edges
      .iter()
      .map(|(&edge, obs)| (edge, obs.fitch_subs().to_vec()))
      .collect(),
  )
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SplitComparison {
  pub shared: usize,
  pub only_first: usize,
  pub only_second: usize,
}

pub struct TreeWithCounts<'a> {
  pub graph: &'a Graph,
  pub names: &'a NodeNames,
  pub edge_counts: &'a BTreeMap<GraphEdgeKey, usize>,
}

pub fn compare_splits(first: &TreeWithCounts<'_>, second: &TreeWithCounts<'_>) -> Result<SplitComparison, Report> {
  let tips = |tree: &TreeWithCounts<'_>| {
    tree
      .graph
      .get_leaves()
      .map(|leaf| tree.names[leaf.key()].to_owned())
      .sorted()
      .collect_vec()
  };
  let first_tips = tips(first);
  if first_tips != tips(second) {
    return make_error!("The two trees have different tip sets");
  }
  let mut ids = TipIds::default();
  for name in &first_tips {
    ids.id(name);
  }
  let first_splits = informative_splits(first, &mut ids, first_tips.len())?;
  let second_splits = informative_splits(second, &mut ids, first_tips.len())?;
  Ok(SplitComparison {
    shared: first_splits.intersection(&second_splits).count(),
    only_first: first_splits.difference(&second_splits).count(),
    only_second: second_splits.difference(&first_splits).count(),
  })
}

fn informative_splits<'a>(
  tree: &TreeWithCounts<'a>,
  ids: &mut TipIds<'a>,
  n_tips: usize,
) -> Result<BTreeSet<Vec<usize>>, Report> {
  let clades = graph_clade_keys(tree.graph, tree.names, ids)?;
  let mut splits = BTreeSet::new();
  for (node, tips) in clades {
    let Some((_, edge)) = tree.graph.node_parent(node)? else {
      continue;
    };
    if tips.len() < 2 || tips.len() + 1 >= n_tips || tree.edge_counts.get(&edge).is_none_or(|&count| count == 0) {
      continue;
    }
    let split = if tips.first() == Some(&0) {
      (0..n_tips).filter(|id| tips.binary_search(id).is_err()).collect_vec()
    } else {
      tips
    };
    splits.insert(split);
  }
  Ok(splits)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MissedSplit {
  pub side: usize,
  pub conflicting: bool,
  pub misplaced: Vec<Displacement>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Displacement {
  pub branches: usize,
  pub substitutions: usize,
}

struct TreeIndex {
  preorder: Vec<GraphNodeKey>,
  parent: BTreeMap<GraphNodeKey, (GraphNodeKey, usize)>,
  clades: BTreeMap<GraphNodeKey, Vec<usize>>,
  tips: BTreeMap<usize, GraphNodeKey>,
  tip_ids: BTreeMap<GraphNodeKey, usize>,
  depth: BTreeMap<GraphNodeKey, Displacement>,
}

impl TreeIndex {
  fn new<'a>(tree: &TreeWithCounts<'a>, ids: &mut TipIds<'a>) -> Result<Self, Report> {
    let clades: BTreeMap<GraphNodeKey, Vec<usize>> =
      graph_clade_keys(tree.graph, tree.names, ids)?.into_iter().collect();
    let schedule = GraphSchedule::new(tree.graph)?;
    let preorder = schedule.preorder().to_vec();
    let mut parent = BTreeMap::new();
    let mut depth = BTreeMap::new();
    for &node in &preorder {
      let here = match tree.graph.node_parent(node)? {
        Some((up, edge)) => {
          let count = tree.edge_counts.get(&edge).copied().unwrap_or(0);
          parent.insert(node, (up, count));
          let above: Displacement = depth[&up];
          Displacement {
            branches: above.branches + 1,
            substitutions: above.substitutions + count,
          }
        },
        None => Displacement {
          branches: 0,
          substitutions: 0,
        },
      };
      depth.insert(node, here);
    }
    let tips: BTreeMap<usize, GraphNodeKey> = clades
      .iter()
      .filter(|&(node, _)| tree.graph.is_leaf(*node))
      .filter_map(|(&node, clade)| clade.first().map(|&id| (id, node)))
      .collect();
    let tip_ids = tips.iter().map(|(&id, &node)| (node, id)).collect();
    Ok(Self {
      preorder,
      parent,
      clades,
      tips,
      tip_ids,
      depth,
    })
  }

  fn informative(&self, node: GraphNodeKey, n_tips: usize) -> bool {
    let size = self.clades[&node].len();
    self.parent.get(&node).is_some_and(|&(_, count)| count > 0) && size >= 2 && size + 2 <= n_tips
  }

  fn distance(&self, mut a: GraphNodeKey, mut b: GraphNodeKey) -> Displacement {
    let (start_a, start_b) = (self.depth[&a], self.depth[&b]);
    while a != b {
      if self.depth[&a].branches >= self.depth[&b].branches {
        a = self.parent.get(&a).map_or(a, |&(up, _)| up);
      } else {
        b = self.parent.get(&b).map_or(b, |&(up, _)| up);
      }
    }
    let ancestor = self.depth[&a];
    Displacement {
      branches: start_a.branches + start_b.branches - 2 * ancestor.branches,
      substitutions: start_a.substitutions + start_b.substitutions - 2 * ancestor.substitutions,
    }
  }
}

pub fn missed_splits(built: &TreeWithCounts<'_>, reference: &TreeWithCounts<'_>) -> Result<Vec<MissedSplit>, Report> {
  let tip_names = built
    .graph
    .get_leaves()
    .map(|leaf| built.names[leaf.key()].to_owned())
    .sorted()
    .collect_vec();
  let n_tips = tip_names.len();
  let mut ids = TipIds::default();
  for name in &tip_names {
    ids.id(name);
  }
  let built_splits = informative_splits(built, &mut ids, n_tips)?;
  let b = TreeIndex::new(built, &mut ids)?;
  let r = TreeIndex::new(reference, &mut ids)?;
  let mut missed = vec![];
  for (&node, clade) in &r.clades {
    if !r.informative(node, n_tips) {
      continue;
    }
    let normalized = if clade.first() == Some(&0) {
      (0..n_tips).filter(|id| clade.binary_search(id).is_err()).collect_vec()
    } else {
      clade.clone()
    };
    if built_splits.contains(&normalized) {
      continue;
    }
    let mut in_side = vec![false; n_tips];
    for &id in clade {
      in_side[id] = true;
    }
    let inside = clade.len();
    let outside = n_tips - inside;
    let mut count: BTreeMap<GraphNodeKey, usize> = BTreeMap::new();
    for &bnode in b.preorder.iter().rev() {
      let own = b.tip_ids.get(&bnode).map_or(0, |&id| usize::from(in_side[id]));
      let total = *count.entry(bnode).or_insert(0) + own;
      count.insert(bnode, total);
      if let Some(&(up, _)) = b.parent.get(&bnode) {
        *count.entry(up).or_insert(0) += total;
      }
    }
    let conflicting = b.preorder.iter().any(|&bnode| {
      if !b.informative(bnode, n_tips) {
        return false;
      }
      let (hit, size) = (count[&bnode], b.clades[&bnode].len());
      hit > 0 && inside > hit && size > hit && outside + hit > size
    });
    let mut misplaced = vec![];
    if conflicting {
      let best = b
        .preorder
        .iter()
        .filter(|bnode| b.parent.contains_key(bnode))
        .map(|&bnode| {
          let (hit, size) = (count[&bnode], b.clades[&bnode].len());
          let same = size + inside - 2 * hit;
          let flipped = 2 * hit + outside - size;
          (same.min(flipped), same <= flipped, bnode)
        })
        .min();
      if let Some((_, same, bnode)) = best {
        let members = &b.clades[&bnode];
        let target_inside = node;
        let target_outside = r.parent.get(&node).map_or(node, |&(up, _)| up);
        for id in 0..n_tips {
          let in_clade = members.binary_search(&id).is_ok();
          let wrong = if same {
            in_clade != in_side[id]
          } else {
            in_clade == in_side[id]
          };
          if wrong {
            let target = if in_side[id] { target_inside } else { target_outside };
            misplaced.push(r.distance(r.tips[&id], target));
          }
        }
      }
    }
    missed.push(MissedSplit {
      side: inside.min(outside),
      conflicting,
      misplaced,
    });
  }
  Ok(missed)
}
