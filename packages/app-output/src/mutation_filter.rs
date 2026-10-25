use eyre::Report;
use std::collections::BTreeMap;
use treetime::ancestral::aa::AaNodeData;
use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::key_map::NodeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
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
    mut edge_mutations: EdgeMap<Vec<Mutation>>,
  ) -> Result<EdgeMap<Vec<Mutation>>, Report> {
    if self.report_unknown {
      return Ok(edge_mutations);
    }
    let schedule = GraphSchedule::new(graph)?;
    schedule.require_tree()?;
    let mut bridge = UnknownBridge::new(self.unknown);
    let mut reported = EdgeMap::new();
    for &key in schedule.preorder() {
      if let &[(parent_key, edge_key)] = schedule.parents(key) {
        let mutations = edge_mutations.remove(edge_key);
        let present = mutations.is_some();
        let mutations = bridge.bridge_edge(
          parent_key,
          key,
          schedule.children(key).len(),
          mutations.unwrap_or_default(),
        )?;
        if present {
          reported.insert(edge_key, mutations);
        }
      }
    }
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
    let schedule = GraphSchedule::new(graph)?;
    schedule.require_tree()?;
    let mut bridge = UnknownBridge::new(self.unknown);
    let mut reported = NodeMap::new();
    for &key in schedule.preorder() {
      let cds_events = node_aa_mutations.remove(key);
      let &[(parent_key, _)] = schedule.parents(key) else {
        if let Some(cds_events) = cds_events {
          reported.insert(key, cds_events);
        }
        continue;
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
      let mutations = bridge.bridge_edge(parent_key, key, schedule.children(key).len(), mutations)?;
      for Mutation { track, event } in mutations {
        if let MutationTrack::AminoAcid(cds) = track {
          by_cds.entry(cds).or_default().push(event);
        }
      }
      if present {
        reported.insert(key, by_cds);
      }
    }
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
    if child_count == 0 {
      return self.bridge_leaf_edge(parent_key, mutations);
    }
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
        reported.extend(bridged_substitution(mutation.track, hidden.remove(&lane), sub)?);
      } else {
        reported.push(mutation);
      }
    }
    if !hidden.is_empty() {
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

  fn bridge_leaf_edge(&mut self, parent_key: GraphNodeKey, mutations: Vec<Mutation>) -> Result<Vec<Mutation>, Report> {
    let parent_states = self.hidden.get(&parent_key).map(|parent| &parent.states);
    let mut reported = Vec::with_capacity(mutations.len());
    for mutation in mutations {
      let MutationEvent::Substitution(sub) = &mutation.event else {
        reported.push(mutation);
        continue;
      };
      if sub.qry() == self.unknown {
        continue;
      }
      if sub.reff() == self.unknown {
        let lane = (mutation.track.clone(), sub.pos());
        let state = parent_states.and_then(|states| states.get(&lane)).copied();
        reported.extend(bridged_substitution(mutation.track, state, sub)?);
      } else {
        reported.push(mutation);
      }
    }
    self.release(parent_key);
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

  fn release(&mut self, parent_key: GraphNodeKey) {
    if let Some(parent) = self.hidden.get_mut(&parent_key) {
      parent.pending_children -= 1;
      if parent.pending_children == 0 {
        self.hidden.remove(&parent_key);
      }
    }
  }
}

struct HiddenStates {
  states: BTreeMap<(MutationTrack, usize), AsciiChar>,
  pending_children: usize,
}

fn bridged_substitution(
  track: MutationTrack,
  hidden: Option<AsciiChar>,
  sub: &Sub,
) -> Result<Option<Mutation>, Report> {
  hidden
    .filter(|&state| state != sub.qry())
    .map(|state| Ok(Mutation::substitution(track, Sub::new(state, sub.pos(), sub.qry())?)))
    .transpose()
}
