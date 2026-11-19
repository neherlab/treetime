#[cfg(test)]
mod tests {
  use crate::optimize::dispatch::guess_initial_branch_lengths;
  use crate::optimize::gather::InitialGuessEvidence;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::seq::indel::{InDel, InDelKind};
  use crate::test_utils::dense_reconstruction_mut;
  use eyre::Report;
  use helpers::setup_dense;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use treetime_primitives::seq::Seq;

  const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  #[test]
  fn test_initial_guess_keep_preserves_zero_lengths_without_indels() -> Result<(), Report> {
    let (graph, reconstruction, mut branch_lengths) = setup_dense(TREE_ZERO_BL)?;

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths)?;
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Keep,
      &mut branch_lengths,
    )?;

    let actual = graph
      .get_edges()
      .map(|edge_ref| branch_lengths[&edge_ref.key()])
      .collect_vec();
    assert_eq!(vec![Some(0.0); 4], actual);
    Ok(())
  }

  #[test]
  fn test_initial_guess_keep_replaces_zero_length_only_on_edge_with_indel() -> Result<(), Report> {
    let (graph, mut reconstruction, mut branch_lengths) = setup_dense(TREE_ZERO_BL)?;

    let edge_key = graph.get_edges().next().expect("tree has edges").key();
    {
      let partition = dense_reconstruction_mut(&mut reconstruction);
      partition
        .edges
        .estimates
        .entry(edge_key)
        .or_default()
        .indels
        .push(InDel {
          range: (4, 7),
          seq: Seq::try_from_str("ACG")?,
          kind: InDelKind::Deletion,
          split_from: None,
        });
    }

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths)?;
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Keep,
      &mut branch_lengths,
    )?;

    let actual = graph
      .get_edges()
      .map(|edge_ref| branch_lengths[&edge_ref.key()])
      .collect_vec();
    let (indel_edge, other_edges) = actual.split_first().expect("tree has edges");
    assert!(
      indel_edge.is_some_and(|bl| bl > 0.0),
      "Keep must replace the zero length of the edge with an indel, got {indel_edge:?}"
    );
    assert_eq!(vec![Some(0.0); 3], other_edges.to_vec());
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::marginal::reconstruction::DenseReconstruction;
    use crate::partition::set::PartitionSet;
    use crate::test_utils::leaf_seq_inputs;
    use eyre::Report;
    use indoc::indoc;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::EdgeMap;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) fn setup_dense(newick: &str) -> Result<(Graph, PartitionSet, EdgeMap<Option<f64>>), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = fasta_read(
        indoc! {r#"
          >A
          AAAACCCCGGGGTTTT
          >B
          CCCCGGGGTTTTAAAA
          >C
          GGGGTTTTAAAACCCC
        "#}
        .as_bytes(),
        &alphabet,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let partition = PartitionMarginalDense::new(alphabet, &graph, &leaf_seq_inputs(&graph, &names, aln))?;
      let reconstruction = PartitionSet::Dense(vec![DenseReconstruction::seeded(
        partition,
        jc69(JC69Params::default())?,
      )]);

      let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

      Ok((graph, reconstruction, branch_lengths))
    }
  }
}
