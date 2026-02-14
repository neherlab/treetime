#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::commands::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::pretty_assert_ulps_eq;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::prelude::*;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::json::{JsonPretty, json_write_str};
  use treetime_io::nwk::nwk_read_str;

  fn assert_dense_rows_normalized(dis: &Array2<f64>, max_ulps: u32) {
    for (row_idx, row) in dis.rows().into_iter().enumerate() {
      let sum: f64 = row.sum();
      assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      assert!(
        sum.is_finite(),
        "Row {row_idx} sum={sum} is not normalized to 1.0 within max_ulps={max_ulps}"
      );
      for (col_idx, &val) in row.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Row {row_idx}, col {col_idx} has non-finite value: {val}"
        );
        assert!(val >= -1e-15, "Row {row_idx}, col {col_idx} has negative value: {val}");
      }
    }
  }

  fn make_nonuniform_gtr(treat_gap_as_unknown: bool) -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })
  }

  fn run_dense_marginal(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
    gtr: GTR,
    treat_gap_as_unknown: bool,
  ) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let log_lh = initialize_marginal(graph, &partitions, aln)?;
    Ok((log_lh, partitions))
  }

  fn run_dense_lh_for_newick(
    newick: &str,
    aln: &[FastaRecord],
    gtr: GTR,
    treat_gap_as_unknown: bool,
  ) -> Result<f64, Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let (log_lh, _) = run_dense_marginal(&graph, aln, gtr, treat_gap_as_unknown)?;
    Ok(log_lh)
  }

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
      &*NUC_ALPHABET,
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
      &*NUC_ALPHABET,
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

    initialize_marginal(&graph, &partitions_marginal_dense, &aln)?;

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

  #[test]
  fn test_marginal_dense_probability_normalization() -> Result<(), Report> {
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
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;

    let (_, partitions) = run_dense_marginal(&graph, &aln, gtr, true)?;
    let partition = partitions[0].read_arc();
    let max_ulps = 4;

    for node_data in partition.nodes.values() {
      if !node_data.profile.dis.is_empty() {
        assert_dense_rows_normalized(&node_data.profile.dis, max_ulps);
      }
    }

    for edge_data in partition.edges.values() {
      if !edge_data.msg_to_child.dis.is_empty() {
        assert_dense_rows_normalized(&edge_data.msg_to_child.dis, max_ulps);
      }
    }

    Ok(())
  }

  #[test]
  fn test_marginal_dense_update_is_idempotent() -> Result<(), Report> {
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
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;

    let (_, partitions) = run_dense_marginal(&graph, &aln, gtr, true)?;

    let log_lh_first = update_marginal(&graph, &partitions)?;
    let log_lh_second = update_marginal(&graph, &partitions)?;

    assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_dense_log_lh_root_invariance_reversible_model() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGTACGTACGT
      >B
      ACGTACGTACGTACGA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#},
      &*NUC_ALPHABET,
    )?;

    let treat_gap_as_unknown = true;
    let gtr1 = make_nonuniform_gtr(treat_gap_as_unknown)?;
    let gtr2 = make_nonuniform_gtr(treat_gap_as_unknown)?;
    let gtr3 = make_nonuniform_gtr(treat_gap_as_unknown)?;

    let tree1 = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let tree2 = "(A:0.1,B:0.2,(C:0.2,D:0.12)CD:0.15)AB:0.01;";
    let tree3 = "((A:0.1,B:0.2)AB:0.15,C:0.2,D:0.12)CD:0.01;";

    let log_lh1 = run_dense_lh_for_newick(tree1, &aln, gtr1, treat_gap_as_unknown)?;
    let log_lh2 = run_dense_lh_for_newick(tree2, &aln, gtr2, treat_gap_as_unknown)?;
    let log_lh3 = run_dense_lh_for_newick(tree3, &aln, gtr3, treat_gap_as_unknown)?;

    let epsilon = 1e-6;
    assert_ulps_eq!(log_lh1, log_lh2, epsilon = epsilon);
    assert_ulps_eq!(log_lh1, log_lh3, epsilon = epsilon);
    assert_ulps_eq!(log_lh2, log_lh3, epsilon = epsilon);

    Ok(())
  }

  #[test]
  fn test_total_likelihood_marginal_dense_all_triplets() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let alphabet = Alphabet::default();
    let mut total_lh = 0.0;

    // Use asymmetric GTR model with non-uniform stationary distribution
    let mu = 1.0;
    let pi = array![0.9, 0.06, 0.02, 0.02];
    let gtr = GTR::new(GTRParams {
      alphabet: alphabet.clone(),
      W: None,
      pi,
      mu,
    })?;

    let graph: GraphAncestral = nwk_read_str("((A:0.6,B:0.3):0.1,C:0.2)root:0.001;")?;
    // Generate all possible triplets (4^3 = 64 combinations)
    let states = ['A', 'C', 'G', 'T'];
    for &state_a in &states {
      for &state_b in &states {
        for &state_c in &states {
          // Create alignment with single position containing this triplet
          let aln = read_many_fasta_str(format!(">A\n{state_a}\n>B\n{state_b}\n>C\n{state_c}\n"), &*NUC_ALPHABET)?;

          let partitions_marginal_dense = [Arc::new(RwLock::new(PartitionMarginalDense {
            index: 0,
            gtr: gtr.clone(),
            alphabet: alphabet.clone(),
            length: get_common_length(&aln)?,
            nodes: btreemap! {},
            edges: btreemap! {},
          }))];

          initialize_marginal(&graph, &partitions_marginal_dense, &aln)?;

          let log_lh = update_marginal(&graph, &partitions_marginal_dense)?;
          total_lh += log_lh.exp();
        }
      }
    }

    // since we test all possible triplets, the total likelihood should be 1
    pretty_assert_ulps_eq!(1.0, total_lh, epsilon = 1e-6);
    Ok(())
  }
}
