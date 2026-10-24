use crate::alphabet::alphabet::Alphabet;
use crate::build::mat::Mat;
use crate::build::placement::{Placement, Target, place_exhaustive, place_exhaustive_where, place_in_region};
use crate::build::query::Query;
use eyre::Report;
use itertools::Itertools;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Debug)]
struct RegionReference {
  node: GraphNodeKey,
  state: Vec<(usize, AsciiChar)>,
  frame: Seq,
  isolation: usize,
}

#[derive(Debug)]
pub(crate) struct Regions {
  references: Vec<RegionReference>,
  owner: BTreeMap<GraphNodeKey, usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RegionPlacement {
  pub placement: Placement,
  pub certified: bool,
  pub local_diffs: usize,
}

impl Regions {
  pub(crate) fn new(mat: &Mat, nodes: &[GraphNodeKey]) -> Result<Self, Report> {
    let owner = voronoi_cells(mat, nodes);
    let references = nodes
      .iter()
      .enumerate()
      .map(|(index, &node)| {
        let state = mat.node_state(node);
        let mut frame = mat.reference().clone();
        for &(pos, value) in &state {
          frame[pos] = value;
        }
        let outside = |target: Target| match target {
          Target::Node(other) => owner.get(&other) != Some(&index),
          Target::Edge(edge) => {
            let (parent, child) = mat.edge_endpoints(edge);
            owner.get(&parent) != Some(&index) && owner.get(&child) != Some(&index)
          },
        };
        let query = Query::from_diffs(state.clone());
        let isolation = place_exhaustive_where(mat, &query, None, None, true, outside).map_or(usize::MAX, |p| p.cost);
        Ok(RegionReference {
          node,
          state,
          frame,
          isolation,
        })
      })
      .collect::<Result<Vec<_>, Report>>()?;
    Ok(Self { references, owner })
  }

  pub(crate) fn len(&self) -> usize {
    self.references.len()
  }

  pub(crate) fn place(
    &self,
    mat: &Mat,
    seq: &Seq,
    query: &Query,
    alphabet: &Alphabet,
    polymorphic: &[usize],
  ) -> Result<RegionPlacement, Report> {
    let nearest = self
      .references
      .iter()
      .enumerate()
      .map(|(index, reference)| (query.cost_against(&reference.state, mat.reference()), index))
      .min();
    let Some((distance, index)) = nearest else {
      let placement = place_exhaustive(mat, query, None, None, true)?;
      return Ok(RegionPlacement {
        placement,
        certified: false,
        local_diffs: query.diffs().len(),
      });
    };
    let reference = &self.references[index];
    let local = Query::from_seq(seq, &reference.frame, alphabet);
    let inside = place_in_region(mat, &local, &reference.frame, reference.node, |node| {
      self.owner.get(&node) == Some(&index)
    })?;
    let slack = distance + query.n_missing_in(polymorphic);
    if inside.cost + slack < reference.isolation {
      return Ok(RegionPlacement {
        placement: inside,
        certified: true,
        local_diffs: local.diffs().len(),
      });
    }
    let exact = place_exhaustive(mat, query, None, Some(inside), true)?;
    Ok(RegionPlacement {
      placement: Placement {
        visited: inside.visited + exact.visited,
        ..exact
      },
      certified: false,
      local_diffs: local.diffs().len(),
    })
  }
}

fn voronoi_cells(mat: &Mat, nodes: &[GraphNodeKey]) -> BTreeMap<GraphNodeKey, usize> {
  let mut owner: BTreeMap<GraphNodeKey, usize> = BTreeMap::new();
  let mut heap: BinaryHeap<Reverse<(usize, usize, GraphNodeKey)>> = nodes
    .iter()
    .enumerate()
    .map(|(index, &node)| Reverse((0, index, node)))
    .collect();
  while let Some(Reverse((distance, index, node))) = heap.pop() {
    if owner.contains_key(&node) {
      continue;
    }
    owner.insert(node, index);
    let neighbors = mat.children(node).chain(mat.parent(node)).collect_vec();
    for (neighbor, edge) in neighbors {
      if !owner.contains_key(&neighbor) {
        heap.push(Reverse((distance + mat.edge_subs(edge).len(), index, neighbor)));
      }
    }
  }
  owner
}
