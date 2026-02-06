#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal_unified::{ancestral_reconstruction_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::nwk::nwk_read_str;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::graph_ancestral::GraphAncestral;
  use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
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
  fn test_ancestral_reconstruction_marginal_sparse() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
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
      &NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      TCGGCGCTGTATTG--
      >AB
      ACATCGCTGTA-TG--
      >CD
      TCGGCGGTGTATTG--
    "#},
      &NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions_marginal_sparse = [Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions_marginal_sparse, &aln)?;

    let log_lh = update_marginal(&graph, &partitions_marginal_sparse)?;

    // generate ancestral reconstruction and test against expectation
    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal(&graph, false, &partitions_marginal_sparse, |node, seq| {
      actual.insert(node.name.clone(), seq.to_string());
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    // test overall likelihood
    pretty_assert_ulps_eq!(-55.55428499726621, log_lh, epsilon = 1e-6);
    Ok(())
  }

  //   #[test]
  //   fn test_root_state() -> Result<(), Report> {
  //     rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //     let aln = read_many_fasta_str(
  //       indoc! {r#"
  //       >A
  //       ACATCGCCNNA--GAC
  //       >B
  //       GCATCCCTGTA-NG--
  //       >C
  //       CCGGCGATGTRTTG--
  //       >D
  //       TCGGCCGTGTRTTG--
  //     "#},
  //       &NUC_ALPHABET,
  //     )?;
  //     let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
  //
  //     let alphabet = Alphabet::default();
  //
  //     let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
  //     let partitions = compress_sequences(&graph, partitions)?;
  //
  //     // use non-trivial GTR with non-uniform stationary distribution
  //     let mu = 1.0;
  //     let pi = array![0.2, 0.3, 0.15, 0.35];
  //     let gtr = GTR::new(GTRParams {
  //       alphabet: Alphabet::default(),
  //       W: None,
  //       pi,
  //       mu,
  //     })?;
  //
  //     let partitions = partitions
  //       .into_iter()
  //       .map(|p| PartitionLikelihood::from_fitch(gtr.clone(), p))
  //       .collect_vec();
  //
  //     let log_lh = run_marginal_sparse(&graph, &partitions)?;
  //     // note that this LH is slightly different from dense or python treetime due to
  //     // different handling of ambiguous characters (value from test_scripts/ancestral_sparse.py)
  //     pretty_assert_ulps_eq!(-56.946298878390444, log_lh, epsilon = 1e-6);
  //
  //     // test variable position distribution at the root (pos 0)
  //     let root = &graph
  //       .get_exactly_one_root()?
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .sparse_partitions[0];
  //     let pos: usize = 0;
  //     let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
  //     pretty_assert_ulps_eq!(&root.profile.variable[&pos].dis, &pos_zero_root, epsilon = 1e-6);
  //
  //     // pull out internal node AB for further testing
  //     let node_ab = &graph
  //       .get_node(GraphNodeKey(1))
  //       .unwrap()
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .sparse_partitions[0];
  //
  //     // test variable position distribution at internal node (pos 0)
  //     let pos: usize = 0;
  //     let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
  //     pretty_assert_ulps_eq!(&node_ab.profile.variable[&pos].dis, &pos_zero_ab, epsilon = 1e-6);
  //
  //     // test variable position distribution at internal node (pos 3)
  //     let dis_ab = array![
  //       0.0013914677323952813,
  //       0.002087201598592933,
  //       0.042827146239885545,
  //       0.9536941844291262
  //     ];
  //     let pos: usize = 3;
  //     pretty_assert_ulps_eq!(&node_ab.profile.variable[&pos].dis, &dis_ab, epsilon = 1e-6);
  //
  //     Ok(())
  //   }
}
