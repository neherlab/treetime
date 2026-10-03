use std::collections::BTreeMap;
use treetime::ancestral::aa::AaNodeData;
use treetime::seq::mutation::{Mutation, MutationEvent};
use treetime_graph::edge::GraphEdgeKey;
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

  pub(crate) fn is_reported(self, event: &MutationEvent) -> bool {
    self.report_unknown
      || match event {
        MutationEvent::Substitution(substitution) => {
          substitution.reff() != self.unknown && substitution.qry() != self.unknown
        },
        MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => true,
      }
  }

  pub(crate) fn reported_events(self, events: Vec<MutationEvent>) -> Vec<MutationEvent> {
    events.into_iter().filter(|event| self.is_reported(event)).collect()
  }

  pub fn reported_edge_mutations(
    &self,
    edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  ) -> BTreeMap<GraphEdgeKey, Vec<Mutation>> {
    edge_mutations
      .into_iter()
      .map(|(edge_key, mutations)| {
        let mutations = mutations
          .into_iter()
          .filter(|mutation| self.is_reported(&mutation.event))
          .collect();
        (edge_key, mutations)
      })
      .collect()
  }

  pub fn reported_aa_node_data(&self, node_data: AaNodeData) -> AaNodeData {
    let node_aa_mutations = node_data
      .node_aa_mutations
      .into_iter()
      .map(|(node_key, cds_mutations)| {
        let cds_mutations = cds_mutations
          .into_iter()
          .map(|(cds, events)| (cds, self.reported_events(events)))
          .collect();
        (node_key, cds_mutations)
      })
      .collect();
    AaNodeData {
      node_aa_mutations,
      ..node_data
    }
  }
}
