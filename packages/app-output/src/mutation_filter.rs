use eyre::Report;
use std::collections::BTreeMap;
use treetime::ancestral::aa::AaNodeData;
use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AsciiChar;

#[derive(Clone, Copy, Debug)]
pub struct UnknownMutationFilter {
  unknown: AsciiChar,
  report_unknown: bool,
}

impl UnknownMutationFilter {
  pub fn new(unknown: AsciiChar, report_unknown: bool) -> Self {
    Self {
      unknown,
      report_unknown,
    }
  }

  pub fn hiding_unknown(unknown: AsciiChar) -> Self {
    Self::new(unknown, false)
  }

  pub fn reported_edge_mutations(
    self,
    graph: &Graph,
    mut edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  ) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
    if self.report_unknown {
      return Ok(edge_mutations);
    }
    let mut bridge = UnknownBridge::new(self.unknown);
    let mut reported = BTreeMap::new();
    graph.iter_depth_first_preorder_forward(|node| {
      if let [(parent_key, edge_key)] = node.parent_keys.as_slice() {
        let mutations = edge_mutations.remove(edge_key);
        let present = mutations.is_some();
        let mutations = bridge.bridge_edge(
          *parent_key,
          node.key,
          node.child_edge_keys.len(),
          mutations.unwrap_or_default(),
        )?;
        if present {
          reported.insert(*edge_key, mutations);
        }
      }
      Ok(())
    })?;
    Ok(reported)
  }

  pub fn reported_aa_node_data(self, graph: &Graph, node_data: AaNodeData) -> Result<AaNodeData, Report> {
    if self.report_unknown {
      return Ok(node_data);
    }
    let AaNodeData {
      reference,
      mut node_aa_mutations,
      root_aa_sequences,
    } = node_data;
    let mut bridge = UnknownBridge::new(self.unknown);
    let mut reported = BTreeMap::new();
    graph.iter_depth_first_preorder_forward(|node| {
      let cds_events = node_aa_mutations.remove(&node.key);
      let [(parent_key, _)] = node.parent_keys.as_slice() else {
        if let Some(cds_events) = cds_events {
          reported.insert(node.key, cds_events);
        }
        return Ok(());
      };
      let present = cds_events.is_some();
      let cds_events = cds_events.unwrap_or_default();
      let mut by_cds: BTreeMap<String, Vec<MutationEvent>> =
        cds_events.keys().map(|cds| (cds.clone(), vec![])).collect();
      let mutations = cds_events
        .into_iter()
        .flat_map(|(cds, events)| {
          events.into_iter().map(move |event| Mutation {
            track: MutationTrack::AminoAcid(cds.clone()),
            event,
          })
        })
        .collect();
      let mutations = bridge.bridge_edge(*parent_key, node.key, node.child_edge_keys.len(), mutations)?;
      for Mutation { track, event } in mutations {
        if let MutationTrack::AminoAcid(cds) = track {
          by_cds.entry(cds).or_default().push(event);
        }
      }
      if present {
        reported.insert(node.key, by_cds);
      }
      Ok(())
    })?;
    Ok(AaNodeData {
      reference,
      node_aa_mutations: reported,
      root_aa_sequences,
    })
  }
}

pub(crate) struct UnknownBridge {
  unknown: AsciiChar,
  hidden: BTreeMap<GraphNodeKey, HiddenStates>,
}

impl UnknownBridge {
  pub(crate) fn new(unknown: AsciiChar) -> Self {
    Self {
      unknown,
      hidden: BTreeMap::new(),
    }
  }

  pub(crate) fn bridge_edge(
    &mut self,
    parent_key: GraphNodeKey,
    node_key: GraphNodeKey,
    child_count: usize,
    mutations: Vec<Mutation>,
  ) -> Result<Vec<Mutation>, Report> {
    let mut hidden = self.inherit(parent_key);
    let mut reported = Vec::with_capacity(mutations.len());
    for mutation in mutations {
      let MutationEvent::Substitution(sub) = &mutation.event else {
        reported.push(mutation);
        continue;
      };
      let lane = (mutation.track.clone(), sub.pos());
      if sub.qry() == self.unknown {
        hidden.insert(lane, sub.reff());
      } else if sub.reff() == self.unknown {
        if let Some(state) = hidden.remove(&lane).filter(|&state| state != sub.qry()) {
          reported.push(Mutation::substitution(
            mutation.track,
            Sub::new(state, sub.pos(), sub.qry())?,
          ));
        }
      } else {
        reported.push(mutation);
      }
    }
    if child_count > 0 && !hidden.is_empty() {
      self.hidden.insert(
        node_key,
        HiddenStates {
          states: hidden,
          pending_children: child_count,
        },
      );
    }
    Ok(reported)
  }

  fn inherit(&mut self, parent_key: GraphNodeKey) -> BTreeMap<(MutationTrack, usize), AsciiChar> {
    let Some(parent) = self.hidden.get_mut(&parent_key) else {
      return BTreeMap::new();
    };
    parent.pending_children -= 1;
    if parent.pending_children == 0 {
      self
        .hidden
        .remove(&parent_key)
        .map(|parent| parent.states)
        .unwrap_or_default()
    } else {
      parent.states.clone()
    }
  }
}

struct HiddenStates {
  states: BTreeMap<(MutationTrack, usize), AsciiChar>,
  pending_children: usize,
}
