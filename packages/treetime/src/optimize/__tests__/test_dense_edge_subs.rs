#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::storage::dense::{DenseNodeState, DenseSeqDistribution, DenseSeqInfo};
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;
  use treetime_primitives::LogLh;

  #[test]
  fn test_dense_edge_subs_no_false_mutation_from_uniform_outgroup() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_ref = &graph.get_edges().collect::<Vec<_>>()[0];
    let edge_key = edge_ref.key();
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    let parent_posterior = array![[0.0, 0.0, 1.0, 0.0]];
    let child_posterior = array![[0.0, 0.0, 1.0, 0.0]];

    let partition = PartitionMarginalDense::new(0, Alphabet::new(AlphabetName::Nuc)?, 1);
    let node_states = btreemap! {
      parent_key => DenseNodeState {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDistribution::new(parent_posterior, LogLh::ZERO),
      },
      child_key => DenseNodeState {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDistribution::new(child_posterior, LogLh::ZERO),
      },
    };

    let subs = partition.edge_subs(&node_states, &graph, edge_key)?;

    assert_eq!(Vec::<Sub>::new(), subs);
    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_detects_real_mutation_hidden_by_edge_messages() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_ref = &graph.get_edges().collect::<Vec<_>>()[0];
    let edge_key = edge_ref.key();
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    let parent_posterior = array![[1.0, 0.0, 0.0, 0.0]];
    let child_posterior = array![[0.0, 1.0, 0.0, 0.0]];

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let parent_state = alphabet.char(0);
    let child_state = alphabet.char(1);

    let partition = PartitionMarginalDense::new(0, alphabet, 1);
    let node_states = btreemap! {
      parent_key => DenseNodeState {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDistribution::new(parent_posterior, LogLh::ZERO),
      },
      child_key => DenseNodeState {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDistribution::new(child_posterior, LogLh::ZERO),
      },
    };

    let subs = partition.edge_subs(&node_states, &graph, edge_key)?;

    let expected = vec![Sub::new(parent_state, 0_usize, child_state)?];
    assert_eq!(expected, subs);
    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_match_reconstructed_branch_differences() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = divergent_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partition = PartitionMarginalDense::new(0, Alphabet::new(AlphabetName::Nuc)?, get_common_length(&aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
    let recon = DenseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let actual_by_edge: BTreeMap<_, _> = graph
      .get_edges()
      .map(|edge_ref| {
        let edge_key = edge_ref.key();
        let subs = recon.partition.edge_subs(&recon.node_states, &graph, edge_key).unwrap();
        (edge_key, subs)
      })
      .collect();

    let expected_by_edge: BTreeMap<_, _> = graph
      .get_edges()
      .map(|edge_ref| {
        let edge_key = edge_ref.key();
        let parent_key = edge_ref.source();
        let child_key = edge_ref.target();
        let expected = helpers::diff_map_states(
          &recon.partition.alphabet,
          &recon.node_states[&parent_key],
          &recon.node_states[&child_key],
        );
        (edge_key, expected)
      })
      .collect();

    assert_eq!(expected_by_edge, actual_by_edge);
    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_excludes_gap_positions_with_posteriors() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2):0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_ref = &graph.get_edges().collect::<Vec<_>>()[0];
    let edge_key = edge_ref.key();
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    #[rustfmt::skip]
    let parent_posterior = array![
      [1.0,  0.0,  0.0,  0.0],
      [0.0,  1.0,  0.0,  0.0],
      [0.0,  0.0,  1.0,  0.0],
      [1.0,  0.0,  0.0,  0.0],
    ];
    #[rustfmt::skip]
    let child_posterior = array![
      [1.0,  0.0,  0.0,  0.0],
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
      [0.0,  1.0,  0.0,  0.0],
    ];

    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(0, alphabet.clone(), 4);
    let node_states = btreemap! {
      parent_key => DenseNodeState {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDistribution::new(parent_posterior, LogLh::ZERO),
      },
      child_key => DenseNodeState {
        seq: DenseSeqInfo { gaps: vec![(1, 3)], non_char: vec![(1, 3)], ..Default::default() },
        profile: DenseSeqDistribution::new(child_posterior, LogLh::ZERO),
      },
    };

    let subs = partition.edge_subs(&node_states, &graph, edge_key)?;

    let expected = vec![Sub::new(alphabet.char(0), 3_usize, alphabet.char(1))?];
    assert_eq!(expected, subs);
    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_is_canonical_filter_present() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = divergent_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partition = PartitionMarginalDense::new(0, Alphabet::new(AlphabetName::Nuc)?, get_common_length(&aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
    let recon = DenseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let subs = recon.partition.edge_subs(&recon.node_states, &graph, edge_key)?;
      for sub in &subs {
        assert!(
          recon.partition.alphabet.is_canonical(sub.reff()),
          "parent state must be canonical: {sub:?}"
        );
        assert!(
          recon.partition.alphabet.is_canonical(sub.qry()),
          "child state must be canonical: {sub:?}"
        );
      }
    }
    Ok(())
  }

  fn divergent_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
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
      &alphabet,
    )
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::partition::storage::dense::DenseNodeState;
    use crate::seq::mutation::Sub;
    use treetime_utils::array::ndarray::argmax_first;

    pub fn diff_map_states(alphabet: &Alphabet, parent_node: &DenseNodeState, child_node: &DenseNodeState) -> Vec<Sub> {
      let parent_profile = &parent_node.profile.dis;
      let child_profile = &child_node.profile.dis;
      let parent_gaps = &parent_node.seq.gaps;
      let child_gaps = &child_node.seq.gaps;

      let mut subs = Vec::new();
      for (pos, (parent_row, child_row)) in parent_profile.rows().into_iter().zip(child_profile.rows()).enumerate() {
        if parent_gaps.iter().any(|&(start, end)| pos >= start && pos < end) {
          continue;
        }
        if child_gaps.iter().any(|&(start, end)| pos >= start && pos < end) {
          continue;
        }

        let parent_state = alphabet.char(argmax_first(&parent_row).unwrap_or(0));
        let child_state = alphabet.char(argmax_first(&child_row).unwrap_or(0));
        if parent_state != child_state && alphabet.is_canonical(parent_state) && alphabet.is_canonical(child_state) {
          subs.push(Sub::new(parent_state, pos, child_state).unwrap());
        }
      }
      subs
    }
  }
}
