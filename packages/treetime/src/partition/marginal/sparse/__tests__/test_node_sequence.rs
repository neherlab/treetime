#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
  use crate::partition::marginal::sparse::reconstruct::LeafDowns;
  use crate::test_utils::find_node_key_by_name;
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::assert_error;

  #[test]
  fn test_sparse_node_sequence_leaf_with_two_parents_errors() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.1)X:0.1,C:0.1)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let mut graph = nwk_parsed.graph;
    let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
    let aln = [("A", "ACGT"), ("B", "ACGT"), ("C", "ACGA")]
      .into_iter()
      .map(|(name, seq)| {
        Ok(AlignmentRecord {
          name: name.to_owned(),
          seq: Seq::try_from_str(seq)?,
        })
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let fitch = create_fitch_partition(&graph, Alphabet::default(), leaf_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let MarginalUpdate { node_states, edges, .. } =
      partition.marginal_update(&jc69(JC69Params::default())?, &graph, &branch_lengths, &node_states)?;
    let x_key = find_node_key_by_name(&graph, &names, "X").expect("node X must exist");
    let c_key = find_node_key_by_name(&graph, &names, "C").expect("leaf C must exist");
    graph.add_edge(x_key, c_key)?;

    let result = partition.node_sequence(
      &graph,
      &node_states,
      LeafDowns::Single {
        forward: &edges.forward,
        node_states: &node_states,
      },
      false,
      c_key,
    );

    assert_error!(
      result,
      format!("Node {c_key} has 2 parents that pass on all sites, but every site has one parent")
    );
    Ok(())
  }
}
