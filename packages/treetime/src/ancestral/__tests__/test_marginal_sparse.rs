#![allow(
  clippy::disallowed_methods,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::seq::alignment::node_seq_inputs;

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::{ancestral_reconstruction, branch_lengths_or_zero};
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::ancestral::sample::SampleMode;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::{GTR, GTRParams};
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::SparseSeqDistribution;
  use crate::pretty_assert_ulps_eq;
  use crate::seq::composition::Composition;
  use crate::seq::mutation::Sub;
  use crate::test_utils::find_node_key_by_name;
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
  use treetime_primitives::{AlignmentRecord, AlphabetLike, Seq};
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  fn assert_sparse_profile_normalized(profile: &SparseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.value().is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh.value()
    );

    for (pos, var_pos) in &profile.variable {
      let sum: f64 = var_pos.dis.sum();
      pretty_assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      assert!(
        sum.is_finite(),
        "Variable position {pos} sum={sum} is not normalized to 1.0 within max_ulps={max_ulps}"
      );
      for (idx, &val) in var_pos.dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Variable position {pos}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Variable position {pos}, index {idx} has negative value: {val}"
        );
      }
    }

    for (char_key, fixed_dis) in &profile.fixed {
      let sum: f64 = fixed_dis.sum();
      pretty_assert_ulps_eq!(sum, 1.0, max_ulps = max_ulps);
      assert!(
        sum.is_finite(),
        "Fixed distribution for char {char_key:?} sum={sum} is not normalized to 1.0 within max_ulps={max_ulps}"
      );
      for (idx, &val) in fixed_dis.iter().enumerate() {
        assert!(
          val.is_finite(),
          "Fixed distribution for char {char_key:?}, index {idx} has non-finite value: {val}"
        );
        assert!(
          val >= -1e-15,
          "Fixed distribution for char {char_key:?}, index {idx} has negative value: {val}"
        );
      }
    }
  }

  fn make_nonuniform_gtr() -> Result<GTR, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    GTR::new(GTRParams {
      n_states,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.35],
    })
  }

  fn run_sparse_marginal(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    gtr: GTR,
  ) -> Result<(f64, SparseReconstruction), Report> {
    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(graph, 0, alphabet, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);
    let (recon, log_lh) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  fn run_sparse_lh_for_newick(newick: &str, aln: &[AlignmentRecord], gtr: GTR) -> Result<f64, Report> {
    let nwk_parsed = nwk_read_str(newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (log_lh, _) = run_sparse_marginal(&graph, &branch_lengths, &names, aln, gtr)?;
    Ok(log_lh)
  }

  #[test]
  fn test_ancestral_reconstruction_marginal_sparse() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

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

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;

    let alphabet = Alphabet::default();

    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);

    let (mut recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();

    let mut actual = BTreeMap::new();
    {
      let SparseReconstruction {
        partition,
        node_states,
        edges,
        ..
      } = &mut recon;
      let mut rng = rand::thread_rng();
      ancestral_reconstruction(&graph, |node| {
        let seq = partition.reconstruct_node_sequence(
          node_states,
          &edges.forward,
          node,
          false,
          false,
          SampleMode::Argmax,
          &mut rng,
        )?;
        actual.insert(names[&node.key].clone(), seq.to_string());
        Some(())
      })?;
    }

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    for name in expected.keys() {
      let node_key = find_node_key_by_name(&graph, &names, name).expect("expected internal node must exist");
      let sequence = &recon.node_states[&node_key].sequence;
      let stored_composition = Composition::with_seq(
        sequence,
        recon.partition.alphabet.chars(),
        recon.partition.alphabet.gap(),
      );
      assert_eq!(stored_composition, recon.partition.obs_nodes[&node_key].composition);
    }

    pretty_assert_ulps_eq!(-55.33813399214274, log_lh, epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_marginal_sparse_probability_normalization() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let gtr = jc69(JC69Params::default())?;

    let (log_lh, recon) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr)?;

    pretty_assert_ulps_eq!(-55.33813399214274, log_lh, epsilon = 1e-6);

    for node_data in recon.node_states.values() {
      assert_sparse_profile_normalized(&node_data.profile, 4);
    }

    for edge_data in recon.edges.forward.values() {
      assert_sparse_profile_normalized(&edge_data.msg_to_child, 4);
    }

    Ok(())
  }

  #[test]
  fn test_marginal_sparse_update_is_idempotent() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let gtr = jc69(JC69Params::default())?;

    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);

    let (recon, log_lh_first) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_first = log_lh_first.value();
    let (recon, log_lh_second) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh_second = log_lh_second.value();

    pretty_assert_ulps_eq!(-55.33813399214274, log_lh_first, epsilon = 1e-6);

    pretty_assert_ulps_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_marginal_sparse_log_lh_root_invariance_reversible_model() -> Result<(), Report> {
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

    let log_lh1 = run_sparse_lh_for_newick(tree1, &aln, gtr1)?;
    let log_lh2 = run_sparse_lh_for_newick(tree2, &aln, gtr2)?;
    let log_lh3 = run_sparse_lh_for_newick(tree3, &aln, gtr3)?;

    pretty_assert_ulps_eq!(log_lh1, log_lh2, epsilon = 1e-6);
    pretty_assert_ulps_eq!(log_lh1, log_lh3, epsilon = 1e-6);
    pretty_assert_ulps_eq!(log_lh2, log_lh3, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_marginal_sparse_posterior_values_python_parity() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let gtr = make_nonuniform_gtr()?;

    let (log_lh, recon) = run_sparse_marginal(&graph, &branch_lengths, &names, &aln, gtr)?;

    pretty_assert_ulps_eq!(-56.76471324493305, log_lh, epsilon = 1e-6);

    let root_key = graph.get_exactly_one_root()?.key();
    let root_profile = &recon.node_states[&root_key].profile;
    let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
    pretty_assert_ulps_eq!(&root_profile.variable[&0].dis, &pos_zero_root, epsilon = 1e-6);

    let ab_key = find_node_key_by_name(&graph, &names, "AB").expect("AB node should exist");
    let ab_profile = &recon.node_states[&ab_key].profile;
    let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
    pretty_assert_ulps_eq!(&ab_profile.variable[&0].dis, &pos_zero_ab, epsilon = 1e-6);

    let dis_ab_pos3 = array![
      0.0013914677323952813,
      0.002087201598592933,
      0.042827146239885545,
      0.9536941844291262
    ];
    pretty_assert_ulps_eq!(&ab_profile.variable[&3].dis, &dis_ab_pos3, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_total_likelihood_marginal_sparse_all_triplets() -> Result<(), Report> {
    let alphabet = Alphabet::default();
    let mut total_lh = 0.0;

    let mu = 1.0;
    let pi = array![0.9, 0.06, 0.02, 0.02];
    let gtr = GTR::new(GTRParams {
      n_states: alphabet.n_canonical(),
      W: None,
      pi,
      mu,
    })?;

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

          let fitch = create_fitch_partition(
            &graph,
            0,
            alphabet.clone(),
            &node_seq_inputs(&graph, &names, aln.clone()),
          )?;
          let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
          let recon = SparseReconstruction::seeded(partition, gtr.clone(), node_states);

          let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
          let log_lh = log_lh.value();
          total_lh += log_lh.exp();
        }
      }
    }

    pretty_assert_ulps_eq!(1.0, total_lh, epsilon = 1e-6);
    Ok(())
  }

  #[test]
  fn test_sparse_edge_subs_match_reconstructed_branch_differences() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let fitch = create_fitch_partition(&graph, 0, Alphabet::default(), &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, make_nonuniform_gtr()?, node_states);
    let (mut recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let actual_by_edge = {
      graph
        .get_edges()
        .map(|edge| {
          let edge_key = edge.key();
          let actual = recon.edge_subs(edge_key)?;
          Ok((edge_key, actual))
        })
        .collect::<Result<BTreeMap<_, _>, Report>>()?
    };

    let mut seqs_by_name = BTreeMap::new();
    {
      let SparseReconstruction {
        partition,
        node_states,
        edges,
        ..
      } = &mut recon;
      let mut rng = rand::thread_rng();
      ancestral_reconstruction(&graph, |node| {
        let seq = partition.reconstruct_node_sequence(
          node_states,
          &edges.forward,
          node,
          true,
          false,
          SampleMode::Argmax,
          &mut rng,
        )?;
        seqs_by_name.insert(names[&node.key].clone().expect("all test nodes should have names"), seq);
        Some(())
      })?;
    }

    let expected_by_edge = helpers::expected_edge_subs_by_edge(&graph, &names, &recon.partition, &seqs_by_name)?;

    assert_eq!(expected_by_edge, actual_by_edge);
    Ok(())
  }

  #[test]
  fn test_marginal_sparse_parallel_pipeline_is_thread_count_deterministic() -> Result<(), Report> {
    let expected = helpers::run_thread_determinism_case(1)?;
    let actual = helpers::run_thread_determinism_case(4)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use super::*;
    use rayon::ThreadPoolBuilder;

    pub fn run_thread_determinism_case(threads: usize) -> Result<(u64, String, String), Report> {
      let newick = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
      let alignment: Vec<AlignmentRecord> = read_many_fasta_str(
        indoc! {r#"
          >A
          ACAACG
          >B
          ACGACG
          >C
          TCGTCG
          >D
          TCGTAG
        "#},
        &*NUC_ALPHABET,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      ThreadPoolBuilder::new().num_threads(threads).build()?.install(|| {
        let nwk_parsed = nwk_read_str(newick)?;
        let names = nwk_parsed.names();
        let graph = nwk_parsed.graph;
        let branch_lengths = nwk_parsed.branch_lengths;
        let (_, recon) = run_sparse_marginal(
          &graph,
          &branch_lengths,
          &names,
          &alignment,
          jc69(JC69Params::default())?,
        )?;
        Ok((
          recon.node_states[&graph.get_exactly_one_root()?.key()]
            .profile
            .log_lh
            .value()
            .to_bits(),
          json_write_str(&recon.node_states, JsonPretty(false))?,
          json_write_str(
            &(&recon.edges.backward, &recon.edges.forward, &recon.edges.estimates),
            JsonPretty(false),
          )?,
        ))
      })
    }

    fn diff_canonical_subs(alphabet: &Alphabet, parent_seq: &Seq, child_seq: &Seq) -> Result<Vec<Sub>, Report> {
      parent_seq
        .iter()
        .zip(child_seq.iter())
        .enumerate()
        .filter(|(_, (parent, child))| parent != child)
        .filter(|(_, (parent, child))| alphabet.is_canonical(**parent) && alphabet.is_canonical(**child))
        .map(|(pos, (parent, child))| Sub::new(*parent, pos, *child))
        .collect()
    }

    pub fn expected_edge_subs_by_edge(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      partition: &PartitionMarginalSparse,
      seqs_by_name: &BTreeMap<String, Seq>,
    ) -> Result<BTreeMap<GraphEdgeKey, Vec<Sub>>, Report> {
      graph
        .get_edges()
        .map(|edge| {
          let edge_key = edge.key();
          let expected = diff_canonical_subs(
            &partition.alphabet,
            get_reconstructed_seq(graph, names, seqs_by_name, edge.source()),
            get_reconstructed_seq(graph, names, seqs_by_name, edge.target()),
          )?;
          Ok((edge_key, expected))
        })
        .collect()
    }

    fn get_reconstructed_seq<'a>(
      _graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      seqs_by_name: &'a BTreeMap<String, Seq>,
      node_key: GraphNodeKey,
    ) -> &'a Seq {
      let name = names[&node_key].as_ref().expect("all test nodes should have names");
      &seqs_by_name[name]
    }
  }
}
