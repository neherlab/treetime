#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
  use crate::seq::alignment::node_seq_inputs;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::assert_error;

  #[test]
  fn test_sparse_node_sequence_leaf_with_two_parents_errors() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.1)X:0.1,C:0.1)root;".as_slice())?;
    let names = nwk_parsed.names();
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
    let fitch = create_fitch_partition(&graph, 0, Alphabet::default(), node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let MarginalUpdate { node_states, edges, .. } =
      partition.marginal_update(&jc69(JC69Params::default())?, &graph, &branch_lengths, &node_states)?;
    let x_key = find_node_key_by_name(&graph, &names, "X").expect("node X must exist");
    let c_key = find_node_key_by_name(&graph, &names, "C").expect("leaf C must exist");
    graph.add_edge(x_key, c_key)?;

    let result = partition.node_sequence(&graph, &node_states, &edges.forward, false, c_key);

    assert_error!(
      result,
      format!(
        "When reconstructing the sequence of node {c_key}: Only trees with exactly one parent per node are currently \
         supported, but node '{c_key}' has 2 parents. This is an internal error. Please report it to developers."
      )
    );
    Ok(())
  }
}
