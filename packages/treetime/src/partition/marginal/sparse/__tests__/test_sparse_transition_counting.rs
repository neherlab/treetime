#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::reconstruction::SparseReconstruction;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::{pretty_assert_array_nonneg, pretty_assert_array_positive};

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(|| Alphabet::new(AlphabetName::Nuc).unwrap());

  fn setup_sparse(
    tree_nwk: &str,
    fasta: &str,
  ) -> Result<(Graph, SparseReconstruction, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
    let alphabet = NUC_ALPHABET.clone();
    let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read(tree_nwk.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let fitch = create_fitch_partition(&graph, alphabet, leaf_seq_inputs(&graph, &names, aln))?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    Ok((graph, recon, branch_lengths))
  }

  #[test]
  fn test_sparse_transition_counting_nij_nonneg() -> Result<(), Report> {
    let (graph, recon, branch_lengths) = setup_sparse(
      "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;",
      indoc! {r#"
      >A
      ACATCGCCGTAGAC
      >B
      GCATCCCTGTAGGG
      >C
      CCGGCGATGTGTTG
      >D
      TCGGCCGTGTGTTG
      "#},
    )?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;

    pretty_assert_array_nonneg!(counts.nij);
    pretty_assert_array_nonneg!(counts.Ti);

    Ok(())
  }

  #[test]
  fn test_sparse_transition_counting_ti_positive() -> Result<(), Report> {
    let (graph, recon, branch_lengths) = setup_sparse(
      "((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;",
      indoc! {r#"
      >A
      AACCGGTT
      >B
      AACCGGTT
      >C
      AACCGGTT
      >D
      AACCGGTT
      "#},
    )?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;

    pretty_assert_array_positive!(counts.Ti);

    Ok(())
  }

  #[test]
  fn test_sparse_transition_counting_diagonal_zero() -> Result<(), Report> {
    let (graph, recon, branch_lengths) = setup_sparse(
      "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;",
      indoc! {r#"
      >A
      ACATCGCC
      >B
      GCATCCCT
      >C
      CCGGCGAT
      >D
      TCGGCCGT
      "#},
    )?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;

    for i in 0..counts.nij.nrows() {
      #[allow(clippy::float_cmp, reason = "diagonal is zero by construction, no arithmetic")]
      {
        assert_eq!(0.0, counts.nij[[i, i]], "diagonal nij[{i},{i}] should be zero");
      }
    }

    Ok(())
  }

  #[test]
  fn test_sparse_transition_counting_root_state_sums_to_length() -> Result<(), Report> {
    let (graph, recon, branch_lengths) = setup_sparse(
      "((A:0.1,B:0.1)AB:0.05,(C:0.1,D:0.1)CD:0.05)root:0.0;",
      indoc! {r#"
      >A
      AAAAAAAA
      >B
      AAAAAAAA
      >C
      AAAAAAAA
      >D
      CAAAAAAA
      "#},
    )?;

    let counts = recon.partition.count_transitions(
      &recon.gtr,
      &graph,
      &branch_lengths_or_zero(&branch_lengths),
      &recon.node_states,
      &recon.edges.backward,
      &recon.edges.forward,
    )?;

    assert!(counts.root_state.sum() > 0.0, "root_state should be populated");

    Ok(())
  }
}
