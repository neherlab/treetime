use crate::build::mat::Mat;
use crate::build::placement::{Target, shared_positions};
use crate::build::query::Query;
use crate::make_internal_report;
use eyre::Report;
use smallvec::{SmallVec, smallvec};
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Insertee {
  Sample(usize),
  Subtree(GraphNodeKey),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Inserted {
  pub attached_to: GraphNodeKey,
  pub created: Option<GraphNodeKey>,
  pub touched: SmallVec<[GraphNodeKey; 3]>,
  pub cost: usize,
}

pub(crate) fn insert(mat: &mut Mat, query: &Query, target: Target, insertee: Insertee) -> Result<Inserted, Report> {
  let (attached_to, created, touched, state) = match target {
    Target::Node(node) => (node, None, smallvec![node], mat.node_state(node)),
    Target::Edge(edge) => {
      let (parent, child) = mat.edge_endpoints(edge);
      let subs = mat.edge_subs(edge);
      let upper = shared_positions(subs, query, mat.reference());
      if upper.is_empty() {
        (parent, None, smallvec![parent], mat.node_state(parent))
      } else if upper.len() == subs.len() && !mat.is_leaf(child) {
        (child, None, smallvec![child], mat.node_state(child))
      } else {
        let mut state = mat.node_state(parent);
        let upper_subs = subs.iter().filter(|sub| upper.binary_search(&sub.pos()).is_ok());
        mat.apply_subs(&mut state, upper_subs);
        let node = mat.split_edge(edge, &upper)?;
        (node, Some(node), smallvec![parent, node, child], state)
      }
    },
  };
  if mat.is_leaf(attached_to) {
    return Err(make_internal_report!(
      "Cannot attach below the sample leaf {attached_to}"
    ));
  }
  let subs = mat.subs_to_query(&state, query)?;
  let cost = subs.len();
  match insertee {
    Insertee::Sample(sample) => {
      mat.add_sample(attached_to, subs, sample)?;
    },
    Insertee::Subtree(subtree) => mat.attach(attached_to, subtree, subs)?,
  }
  Ok(Inserted {
    attached_to,
    created,
    touched,
    cost,
  })
}
