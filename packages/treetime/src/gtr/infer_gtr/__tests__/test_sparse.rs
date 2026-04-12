//! Tests for sparse GTR inference.

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::ancestral_reconstruction_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::infer_gtr::common::{InferGtrOptions, MutationCounts, infer_gtr_impl};
  use crate::gtr::infer_gtr::sparse::get_mutation_counts_sparse;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use ndarray::{Array1, Array2};
  use parking_lot::RwLock;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  fn expected_mutation_counts_from_reconstructed_sequences(
    graph: &GraphAncestral,
    partition: &Arc<RwLock<PartitionMarginalSparse>>,
  ) -> Result<MutationCounts, Report> {
    let alphabet = &partition.read_arc().alphabet;
    let mut reconstructed_by_name = BTreeMap::new();
    ancestral_reconstruction_marginal(graph, true, std::slice::from_ref(partition), |node, seq| {
      reconstructed_by_name.insert(
        node.name.clone().expect("all test nodes should have names"),
        seq.clone(),
      );
      Ok(())
    })?;

    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    let root_seq = get_reconstructed_seq(graph, &reconstructed_by_name, root_key);
    let root_state = Array1::from_iter(
      alphabet
        .canonical()
        .map(|nuc| root_seq.iter().filter(|&&state| state == nuc).count() as f64),
    );

    let n_states = alphabet.n_canonical();
    let mut nij = Array2::zeros((n_states, n_states));
    let mut Ti = Array1::zeros(n_states);

    for edge_ref in graph.get_edges() {
      let edge = edge_ref.read_arc();
      let branch_length = edge.payload().read_arc().branch_length().unwrap_or(0.0);
      let parent_seq = get_reconstructed_seq(graph, &reconstructed_by_name, edge.source());
      let child_seq = get_reconstructed_seq(graph, &reconstructed_by_name, edge.target());

      for (i, nuc) in alphabet.canonical().enumerate() {
        Ti[i] += branch_length * child_seq.iter().filter(|&&state| state == nuc).count() as f64;
      }

      for (&parent_state, &child_state) in parent_seq.iter().zip(child_seq.iter()) {
        if parent_state == child_state {
          continue;
        }
        if !alphabet.is_canonical(parent_state) || !alphabet.is_canonical(child_state) {
          continue;
        }
        let i = alphabet.index(child_state)?;
        let j = alphabet.index(parent_state)?;
        nij[[i, j]] += 1.0;
        Ti[i] -= 0.5 * branch_length;
        Ti[j] += 0.5 * branch_length;
      }
    }

    Ok(MutationCounts { nij, Ti, root_state })
  }

  fn get_reconstructed_seq<'a>(
    graph: &GraphAncestral,
    reconstructed_by_name: &'a BTreeMap<String, Seq>,
    node_key: treetime_graph::node::GraphNodeKey,
  ) -> &'a Seq {
    let node = graph.get_node(node_key).expect("node should exist");
    let node = node.read_arc();
    let payload = node.payload().read_arc();
    let name = payload.name.as_ref().expect("all test nodes should have names");
    &reconstructed_by_name[name]
  }

  #[test]
  fn test_get_mutation_counts_sparse() -> Result<(), Report> {
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
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&partition), &aln)?;
    // GTR inference runs before marginal in the production sparse flow (see
    // `ancestral/run.rs`), so `node.seq.composition` holds the Fitch-resolved
    // composition at inference time. The oracle reconstructs sequences without
    // running marginal, which produces the same Fitch-resolved sequences and
    // therefore matches.

    let counts_actual = get_mutation_counts_sparse(&graph, &partition)?;
    let counts_expected = expected_mutation_counts_from_reconstructed_sequences(&graph, &partition)?;

    pretty_assert_ulps_eq!(counts_expected.nij, counts_actual.nij, epsilon = 1e-9);
    pretty_assert_ulps_eq!(counts_expected.Ti, counts_actual.Ti, epsilon = 1e-7);
    pretty_assert_ulps_eq!(counts_expected.root_state, counts_actual.root_state, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_infer_gtr_sparse() -> Result<(), Report> {
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
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    compress_sequences(&graph, std::slice::from_ref(&partition), &aln)?;
    // GTR inference runs before marginal in production; see test above.

    let counts_actual = get_mutation_counts_sparse(&graph, &partition)?;
    let counts_expected = expected_mutation_counts_from_reconstructed_sequences(&graph, &partition)?;
    let actual = infer_gtr_impl(
      &counts_actual,
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;
    let expected = infer_gtr_impl(
      &counts_expected,
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    pretty_assert_ulps_eq!(&expected.W, &actual.W, epsilon = 1e-7);
    pretty_assert_ulps_eq!(&expected.pi, &actual.pi, epsilon = 1e-7);
    pretty_assert_ulps_eq!(expected.mu, actual.mu, epsilon = 1e-7);

    Ok(())
  }
}
