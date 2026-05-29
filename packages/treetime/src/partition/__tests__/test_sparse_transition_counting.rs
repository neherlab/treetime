#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::update_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::traits::TransitionCounting;
  use crate::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::{pretty_assert_array_nonneg, pretty_assert_array_positive};

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  }

  fn setup_sparse(
    tree_nwk: &str,
    fasta: &str,
  ) -> Result<
    (
      GraphAncestral,
      Arc<RwLock<crate::partition::marginal_sparse::PartitionMarginalSparse>>,
    ),
    Report,
  > {
    let alphabet = NUC_ALPHABET.clone();
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let graph: GraphAncestral = nwk_read_str(tree_nwk)?;
    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let partition = Arc::new(RwLock::new(fitch.into_marginal_sparse(gtr, &graph)?));
    update_marginal(&graph, std::slice::from_ref(&partition))?;
    Ok((graph, partition))
  }

  #[test]
  fn test_sparse_transition_counting_nij_nonneg() -> Result<(), Report> {
    let (graph, partition) = setup_sparse(
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

    let counts = partition.read_arc().count_transitions(&graph)?;

    pretty_assert_array_nonneg!(counts.nij);
    pretty_assert_array_nonneg!(counts.Ti);

    Ok(())
  }

  #[test]
  fn test_sparse_transition_counting_ti_positive() -> Result<(), Report> {
    let (graph, partition) = setup_sparse(
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

    let counts = partition.read_arc().count_transitions(&graph)?;

    pretty_assert_array_positive!(counts.Ti);

    Ok(())
  }

  #[test]
  fn test_sparse_transition_counting_diagonal_zero() -> Result<(), Report> {
    let (graph, partition) = setup_sparse(
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

    let counts = partition.read_arc().count_transitions(&graph)?;

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
    let (graph, partition) = setup_sparse(
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

    let counts = partition.read_arc().count_transitions(&graph)?;

    assert!(counts.root_state.sum() > 0.0, "root_state should be populated");

    Ok(())
  }
}
