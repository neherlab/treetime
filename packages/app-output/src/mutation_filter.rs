use std::collections::BTreeMap;
use treetime::ancestral::aa::AaNodeData;
use treetime::seq::mutation::{Mutation, MutationEvent};
use treetime_graph::edge::GraphEdgeKey;
use treetime_primitives::AsciiChar;

#[derive(Clone, Copy, Debug)]
pub struct AmbiguousMutationFilter {
  ambiguous: AsciiChar,
  report_ambiguous: bool,
}

impl AmbiguousMutationFilter {
  pub fn new(ambiguous: AsciiChar, report_ambiguous: bool) -> Self {
    Self {
      ambiguous,
      report_ambiguous,
    }
  }

  pub fn is_reported(&self, event: &MutationEvent) -> bool {
    self.report_ambiguous
      || match event {
        MutationEvent::Substitution(substitution) => {
          substitution.reff() != self.ambiguous && substitution.qry() != self.ambiguous
        },
        MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => true,
      }
  }

  pub fn reported_events(&self, events: Vec<MutationEvent>) -> Vec<MutationEvent> {
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
