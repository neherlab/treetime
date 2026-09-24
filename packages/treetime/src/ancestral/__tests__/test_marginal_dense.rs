#![allow(
  clippy::disallowed_methods,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::ancestral::marginal::{ancestral_reconstruction, branch_lengths_or_zero};
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::ancestral::sample::SampleMode;
  use crate::ancestral::tip_states::TipStates;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use treetime_graph::graph::Graph;

  use ndarray::prelude::*;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  fn assert_dense_rows_normalized(dis: &Array2<f64>, max_ulps: u32) {
    for (row_idx, row) in dis.rows().into_iter().enumerate() {
      let sum: f64 = row.sum();
      pretty_assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
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

  fn make_nonuniform_gtr() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    GTR::builder()
      .n_states(n_states)
      .mu(1.0)
      .pi(array![0.2, 0.3, 0.15, 0.35])
      .build()
  }

  fn run_dense_marginal(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    gtr: GTR,
  ) -> Result<(f64, DenseReconstruction), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(graph, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, log_lh) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  fn run_dense_lh_for_newick(newick: &str, aln: &[AlignmentRecord], gtr: GTR) -> Result<f64, Report> {
    let nwk_parsed = nwk_read_str(newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (log_lh, _) = run_dense_marginal(&graph, &branch_lengths, &names, aln, gtr)?;
    Ok(log_lh)
  }

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  static TREE_7_TAXON: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  static ALN_7_TAXON: LazyLock<Vec<AlignmentRecord>> = LazyLock::new(|| {
    read_many_fasta_str(
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
    )
    .unwrap()
    .into_iter()
    .map(AlignmentRecord::from)
    .collect()
  });

  #[test]
  fn test_ancestral_reconstruction_marginal_dense() -> Result<(), Report> {
    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      TCGGCGCTGTATTG--
      >AB
      ACATCGCTGTA--G--
      >CD
      TCGGCGGTGTATTG--
    "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let nwk_parsed = nwk_read_str(TREE_7_TAXON)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (_, mut recon) = run_dense_marginal(&graph, &branch_lengths, &names, &ALN_7_TAXON, gtr)?;

    let mut actual = BTreeMap::new();
    {
      let DenseReconstruction {
        partition, node_states, ..
      } = &mut recon;
      let mut rng = rand::thread_rng();
      ancestral_reconstruction(&graph, |node| {
        let Some(seq) = partition.reconstruct_node_sequence(
          node_states,
          node,
          TipStates {
            include_leaves: false,
            impute: false,
          },
          SampleMode::Argmax,
          &mut rng,
        ) else {
          return Ok(false);
        };
        actual.insert(names[&node.key].clone(), seq.to_string());
        Ok(true)
      })?;
    }

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_marginal_dense_probability_normalization() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_7_TAXON)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh, recon) = run_dense_marginal(&graph, &branch_lengths, &names, &ALN_7_TAXON, gtr)?;

    pretty_assert_ulps_eq!(-57.712498930787206, log_lh, epsilon = 1e-6);

    let max_ulps = 4;

    for node_data in recon.node_states.values() {
      if !node_data.profile.dis.is_empty() {
        assert_dense_rows_normalized(&node_data.profile.dis, max_ulps);
      }
    }

    for edge_data in recon.edges.forward.values() {
      if !edge_data.msg_to_child.dis.is_empty() {
        assert_dense_rows_normalized(&edge_data.msg_to_child.dis, max_ulps);
      }
    }

    Ok(())
  }

  #[test]
  fn test_marginal_dense_update_is_idempotent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_7_TAXON)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;

    let (log_lh_init, recon) = run_dense_marginal(&graph, &branch_lengths, &names, &ALN_7_TAXON, gtr)?;

    let (recon, log_lh_first) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_first = log_lh_first.value();
    let (recon, log_lh_second) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_second = log_lh_second.value();

    pretty_assert_ulps_eq!(log_lh_init, log_lh_first, epsilon = 1e-10);
    pretty_assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_dense_log_lh_root_invariance_reversible_model() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let gtr1 = make_nonuniform_gtr()?;
    let gtr2 = make_nonuniform_gtr()?;
    let gtr3 = make_nonuniform_gtr()?;

    let tree1 = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let tree2 = "(A:0.1,B:0.2,(C:0.2,D:0.12)CD:0.15)AB:0.01;";
    let tree3 = "((A:0.1,B:0.2)AB:0.15,C:0.2,D:0.12)CD:0.01;";

    let log_lh1 = run_dense_lh_for_newick(tree1, &aln, gtr1)?;
    let log_lh2 = run_dense_lh_for_newick(tree2, &aln, gtr2)?;
    let log_lh3 = run_dense_lh_for_newick(tree3, &aln, gtr3)?;

    pretty_assert_ulps_eq!(log_lh1, log_lh2, epsilon = 1e-6);
    pretty_assert_ulps_eq!(log_lh1, log_lh3, epsilon = 1e-6);
    pretty_assert_ulps_eq!(log_lh2, log_lh3, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_total_likelihood_marginal_dense_all_triplets() -> Result<(), Report> {
    let alphabet = Alphabet::default();
    let mut total_lh = 0.0;

    let mu = 1.0;
    let pi = array![0.9, 0.06, 0.02, 0.02];
    let gtr = GTR::builder().n_states(alphabet.n_canonical()).pi(pi).mu(mu).build()?;

    let nwk_parsed = nwk_read_str("((A:0.6,B:0.3):0.1,C:0.2)root:0.001;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let states = ['A', 'C', 'G', 'T'];
    for &state_a in &states {
      for &state_b in &states {
        for &state_c in &states {
          let aln: Vec<AlignmentRecord> =
            read_many_fasta_str(format!(">A\n{state_a}\n>B\n{state_b}\n>C\n{state_c}\n"), &*NUC_ALPHABET)?
              .into_iter()
              .map(AlignmentRecord::from)
              .collect();
          let (log_lh, _) = run_dense_marginal(&graph, &branch_lengths, &names, &aln, gtr.clone())?;
          total_lh += log_lh.exp();
        }
      }
    }

    pretty_assert_ulps_eq!(1.0, total_lh, max_ulps = 4);
    Ok(())
  }
}
