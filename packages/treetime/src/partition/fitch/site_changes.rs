use crate::alphabet::alphabet::Alphabet;
use crate::homoplasy::classify::MutationClass;
use crate::homoplasy::site_branches::{SiteBranches, sites_by_branch_count};
use crate::partition::fitch::passes::create_fitch_partition;
use crate::seq::alignment::NodeSeqInput;
use crate::seq::mutation::MutationEvent;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;

pub fn fitch_site_change_counts(
  graph: &Graph,
  alphabet: Alphabet,
  node_inputs: NodeMap<NodeSeqInput>,
) -> Result<Vec<usize>, Report> {
  let partition = create_fitch_partition(graph, alphabet, node_inputs)?;
  let sites = sites_by_branch_count(
    partition
      .edges
      .values()
      .map(|edge| edge.fitch_subs().iter().cloned().map(MutationEvent::Substitution)),
    &partition.alphabet,
    MutationClass::Substitution,
  );
  Ok(sites.into_iter().fold(
    vec![0; partition.length],
    |mut counts, SiteBranches { position, branches }| {
      counts[position] = branches;
      counts
    },
  ))
}
