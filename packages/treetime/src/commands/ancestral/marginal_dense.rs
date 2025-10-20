#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::commands::ancestral::marginal_unified::{ancestral_reconstruction_marginal, run_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::nwk::nwk_read_str;
  use crate::representation::graph_ancestral::GraphAncestral;
  use crate::representation::partition_marginal_dense::PartitionMarginalDense;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_io::json::{JsonPretty, json_write_str};

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_ancestral_reconstruction_marginal_dense() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      TCAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >CD
      CCGGCCATGTATTG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      TCGGCGCTGTATTGAC
      >AB
      ACATCGCTGTA-TGAC
      >CD
      TCGGCGGTGTATTG--
    "#},
      &NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let treat_gap_as_unknown = true;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;

    let partitions_marginal_dense = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    run_marginal(&graph, &partitions_marginal_dense, Some(&aln))?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal(&graph, false, &partitions_marginal_dense, |node, seq| {
      actual.insert(node.name.clone(), seq.to_string());
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  //   #[test]
  //   fn test_root_state() -> Result<(), Report> {
  //     rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //     let aln = read_many_fasta_str(
  //       indoc! {r#"
  //       >root
  //       ACAGCCATGTATTG--
  //       >AB
  //       ACATCCCTGTA-TG--
  //       >A
  //       ACATCGCCNNA--GAC
  //       >B
  //       GCATCCCTGTA-NG--
  //       >CD
  //       CCGGCCATGTATTG--
  //       >C
  //       CCGGCGATGTRTTG--
  //       >D
  //       TCGGCCGTGTRTTG--
  //     "#},
  //       &NUC_ALPHABET,
  //     )?;
  //     let graph: DenseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
  //
  //     let treat_gap_as_unknown = true;
  //     let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
  //
  //     // use non-trivial GTR with non-uniform stationary distribution (tests correct use of transposed matrices)
  //     let mu = 1.0;
  //     let pi = array![0.2, 0.3, 0.15, 0.35];
  //     let gtr = GTR::new(GTRParams {
  //       alphabet: Alphabet::default(),
  //       W: None,
  //       pi,
  //       mu,
  //     })?;
  //
  //     let partitions = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];
  //
  //     let log_lh = run_marginal_dense(&graph, partitions, true)?;
  //     // from test_scripts/ancestral_dense.py
  //     pretty_assert_ulps_eq!(-59.20297892181229, log_lh, epsilon = 1e-6);
  //
  //     // test variable position distribution at the root for position 0 (from test_scripts/ancestral_dense.py)
  //     let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
  //     let root = &graph
  //       .get_exactly_one_root()?
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .dense_partitions[0];
  //     let pos: usize = 0;
  //     pretty_assert_ulps_eq!(root.profile.dis.slice(s![pos, 0..4]), &pos_zero_root, epsilon = 1e-6);
  //
  //     // pull out internal node AB for testing
  //     let node_ab = &graph
  //       .get_node(GraphNodeKey(1))
  //       .unwrap()
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .dense_partitions[0];
  //
  //     // test variable position distribution at internal node (from test_scripts/ancestral_dense.py)
  //     let pos: usize = 0;
  //     let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
  //     pretty_assert_ulps_eq!(node_ab.profile.dis.slice(s![pos, 0..4]), &pos_zero_ab, epsilon = 1e-6);
  //
  //     // test variable position distribution at internal node (obtained from python treetime)
  //     let dis_ab = array![
  //       0.0013914677323952813,
  //       0.002087201598592933,
  //       0.042827146239885545,
  //       0.9536941844291262
  //     ];
  //     let pos: usize = 3;
  //     pretty_assert_ulps_eq!(node_ab.profile.dis.slice(s![pos, 0..4]), &dis_ab, epsilon = 1e-6);
  //
  //     // test whether the log likelihood is the same regardless of the root (here for node AB)
  //     let log_lh_ab = node_ab.profile.log_lh;
  //     pretty_assert_ulps_eq!(log_lh_ab, log_lh, epsilon = 1e-8);
  //
  //     Ok(())
  //   }
}
