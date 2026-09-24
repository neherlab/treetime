#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::timetree::partition::PartitionTimetree;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  #[test]
  fn test_partition_timetree_edge_sub_count_unavailable_before_inference() -> Result<(), Report> {
    let (graph, seeded) = helpers::seeded_sparse()?;
    let partition = PartitionTimetree::Sparse(seeded);
    for edge in graph.get_edges() {
      assert_eq!(None, partition.edge_sub_count(&graph, edge.key())?);
    }
    Ok(())
  }

  #[test]
  fn test_partition_timetree_edge_sub_count_matches_edge_subs_after_inference() -> Result<(), Report> {
    let (graph, seeded) = helpers::seeded_sparse()?;
    let (updated, _) = seeded.marginal_update(&graph, &branch_lengths_or_zero(&helpers::branch_lengths()?))?;
    let partition = PartitionTimetree::Sparse(updated);
    for edge in graph.get_edges() {
      let expected = partition.edge_subs(&graph, edge.key())?.len();
      assert_eq!(Some(expected), partition.edge_sub_count(&graph, edge.key())?);
    }
    Ok(())
  }

  mod helpers {
    use super::*;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;

    const TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    pub(super) fn seeded_sparse() -> Result<(Graph, SparseReconstruction), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = read_many_fasta_str(
        indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
        &alphabet,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let nwk_parsed = nwk_read_str(TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let gtr = jc69(JC69Params::default())?;
      Ok((graph, SparseReconstruction::seeded(partition, gtr, node_states)))
    }

    pub(super) fn branch_lengths() -> Result<BTreeMap<GraphEdgeKey, Option<f64>>, Report> {
      Ok(nwk_read_str(TREE)?.branch_lengths)
    }
  }
}
