#[cfg(test)]
mod tests {
  use crate::ancestral::plan::ReconstructionMethod;
  use crate::gtr::get_gtr::{F81Params, GtrModel, JC69Params};
  use crate::partition::create::Representation;
  use eyre::Report;
  use helpers::{Reconstructed, masked_clade_sequences, reconstruct};
  use maplit::btreemap;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::{o, vec_of_owned};

  const ORACLE_TREE: &str = "(((A:0.1,B:0.1)X:0.1,C:0.1)Y:0.05,(D:0.1,E:0.2)Z:0.3)R;";

  const ORACLE_ALIGNMENT: &str = ">A\nNNNNACGT\n>B\nNNNGACGT\n>C\nTNNTACGA\n>D\nCNTCACGT\n>E\nCNCCACTT\n";

  fn oracle_model() -> GtrModel {
    GtrModel::F81(F81Params {
      pi: Some(array![0.1, 0.2, 0.4, 0.2]),
      ..F81Params::default()
    })
  }

  #[rstest]
  #[trace]
  fn test_internal_nodes_without_descendant_data_take_the_states_of_v0(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[values(false, true)] impute: bool,
  ) -> Result<(), Report> {
    let Reconstructed { sequences, .. } = reconstruct(
      ORACLE_TREE,
      ORACLE_ALIGNMENT,
      ReconstructionMethod::Marginal(representation),
      &oracle_model(),
      impute,
      true,
    )?;
    let internal = sequences
      .into_iter()
      .filter(|(name, _)| ["R", "X", "Y", "Z"].contains(&name.as_str()))
      .collect();
    assert_eq!(
      btreemap! {
        o!("R") => o!("TGTTACGT"),
        o!("X") => o!("TGTGACGT"),
        o!("Y") => o!("TGTTACGT"),
        o!("Z") => o!("CGTCACGT"),
      },
      internal,
      "v0 TreeAnc.infer_ancestral_sequences(marginal=True) with GTR.custom(pi=[0.1, 0.2, 0.4, 0.2, 0.1], alphabet='nuc')"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::equilibrium_a(array![0.4, 0.2, 0.2, 0.2], "AAAA")]
  #[case::equilibrium_t(array![0.1, 0.2, 0.3, 0.4], "TTTT")]
  #[trace]
  fn test_column_without_data_in_any_tip_takes_the_equilibrium_state(
    #[case] pi: Array1<f64>,
    #[case] expected: &str,
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
  ) -> Result<(), Report> {
    let model = GtrModel::F81(F81Params { pi: Some(pi), ..F81Params::default() });
    let Reconstructed { sequences, .. } = reconstruct(
      "((L1:0.1,L2:0.15)X:0.2,L3:0.3)root;",
      ">L1\nNNNN\n>L2\nNNNN\n>L3\nNNNN\n",
      ReconstructionMethod::Marginal(representation),
      &model,
      false,
      true,
    )?;
    assert_eq!((expected, expected), (sequences["root"].as_str(), sequences["X"].as_str()));
    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_equal_equilibrium_frequencies_resolve_to_the_first_state(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
  ) -> Result<(), Report> {
    let Reconstructed { sequences, .. } = reconstruct(
      "((L1:0.1,L2:0.15)X:0.2,L3:0.3)root;",
      ">L1\nACNNGT\n>L2\nAGNNGT\n>L3\nACNNGA\n",
      ReconstructionMethod::Marginal(representation),
      &GtrModel::JC69(JC69Params::default()),
      false,
      true,
    )?;
    assert_eq!(
      ("ACAAGT", "ACAAGT"),
      (sequences["root"].as_str(), sequences["X"].as_str())
    );
    Ok(())
  }

  #[rstest]
  #[case::dense(ReconstructionMethod::Marginal(Representation::Dense))]
  #[case::sparse(ReconstructionMethod::Marginal(Representation::Sparse))]
  #[trace]
  fn test_changes_to_unknown_states_appear_on_each_tip(
    #[case] method: ReconstructionMethod,
    #[values(false, true)] report_unknown: bool,
  ) -> Result<(), Report> {
    let Reconstructed { mutations, sequences } = reconstruct(
      "(((A:0.1,B:0.1)AB:0.1,C:0.1)ABC:0.05,D:0.3)root;",
      ">A\nNNGAAATTTCCCGGG\n>B\nNNNAAATTTCCCGGG\n>C\nATGAAATTTCCCGNN\n>D\nATGAAATTTCCAGNN\n",
      method,
      &GtrModel::JC69(JC69Params::default()),
      false,
      report_unknown,
    )?;
    let expected = if report_unknown {
      btreemap! {
        o!("A") => vec_of_owned!["A1N", "T2N"],
        o!("B") => vec_of_owned!["A1N", "T2N", "G3N"],
        o!("C") => vec_of_owned!["G14N", "G15N"],
        o!("D") => vec_of_owned!["C12A", "G14N", "G15N"],
        o!("AB") => vec![],
        o!("ABC") => vec![],
      }
    } else {
      btreemap! {
        o!("A") => vec![],
        o!("B") => vec![],
        o!("C") => vec![],
        o!("D") => vec_of_owned!["C12A"],
        o!("AB") => vec![],
        o!("ABC") => vec![],
      }
    };
    assert_eq!(expected, mutations);
    assert_eq!(
      ("ATGAAATTTCCCGGG", "ATGAAATTTCCCGGG", "ATGAAATTTCCCGGG"),
      (
        sequences["root"].as_str(),
        sequences["ABC"].as_str(),
        sequences["AB"].as_str()
      ),
      "v0 treetime ancestral --gtr JC69"
    );
    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_dense_and_sparse_agree_on_columns_without_data_in_any_tip(
    #[values(oracle_model(), GtrModel::JC69(JC69Params::default()))] model: GtrModel,
    #[values(false, true)] impute: bool,
    #[values(false, true)] report_unknown: bool,
  ) -> Result<(), Report> {
    let alignment = ">A\nACGTN--TN\n>B\nACGAN--TN\n>C\nTCGTN-NAN\n>D\nCCTCN-NTN\n>E\nCGCCN-NTN\n";
    let reconstruct_with = |representation| {
      reconstruct(
        ORACLE_TREE,
        alignment,
        ReconstructionMethod::Marginal(representation),
        &model,
        impute,
        report_unknown,
      )
    };
    let dense = reconstruct_with(Representation::Dense)?;
    let sparse = reconstruct_with(Representation::Sparse)?;
    assert_eq!(
      (&dense.sequences, &dense.mutations),
      (&sparse.sequences, &sparse.mutations),
      "dense marginal reconstruction, which v0 cross-checks"
    );
    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_dense_and_sparse_agree_on_a_clade_without_data(#[values(false, true)] impute: bool) -> Result<(), Report> {
    let dense = masked_clade_sequences(Representation::Dense, impute)?;
    let sparse = masked_clade_sequences(Representation::Sparse, impute)?;
    assert_eq!(dense.masked_internal, sparse.masked_internal);
    assert_eq!(dense.sequences, sparse.sequences);
    assert_eq!(dense.mutations, sparse.mutations);
    assert!(
      dense.masked_internal.values().all(|seq| !seq.contains('N')),
      "internal nodes of the masked clade carry no unknown state: {:?}",
      dense.masked_internal
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::plan::{ReconstructionMethod, ReconstructionOptions, reconstruct_partition};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModel;
    use crate::partition::alignment_sequences::PartitionSequences;
    use crate::partition::create::Representation;
    use crate::partition::inferred_rate::InferredRate;
    use crate::partition::marginal::sample::SampleMode;
    use crate::progress::NoopProgress;
    use crate::seq::mutation::{Mutation, MutationEvent, Sub};
    use crate::seq::sink::{SeqStream, SeqTrack};
    use crate::test_utils::{RecordingSeqSink, TestDataset, fasta_xz_records, leaf_seq_inputs};
    use eyre::Report;
    use itertools::Itertools;
    use pretty_assertions::assert_eq;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::NodeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};
    use treetime_utils::make_report;
    use treetime_utils::sync::random::get_random_number_generator;

    pub(super) struct Reconstructed {
      pub mutations: BTreeMap<String, Vec<String>>,
      pub sequences: BTreeMap<String, String>,
    }

    pub(super) fn reconstruct(
      newick: &str,
      fasta: &str,
      method: ReconstructionMethod,
      model: &GtrModel,
      impute: bool,
      report_unknown: bool,
    ) -> Result<Reconstructed, Report> {
      let alphabet = Alphabet::default();
      let sequences = fasta_read(fasta.as_bytes(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect_vec();
      reconstruct_records(newick, sequences, method, model, impute, report_unknown)
    }

    fn reconstruct_records(
      newick: &str,
      sequences: Vec<AlignmentRecord>,
      method: ReconstructionMethod,
      model: &GtrModel,
      impute: bool,
      report_unknown: bool,
    ) -> Result<Reconstructed, Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read(newick.as_bytes())?;
      let names = parse.names.clone();
      let graph = &parse.graph;
      let partition = reconstruct_partition(
        graph,
        method,
        model,
        InferredRate::Unit,
        alphabet.clone(),
        leaf_seq_inputs(graph, &names, sequences),
        &branch_lengths_or_zero(&parse.branch_lengths),
        &ReconstructionOptions::new(impute, SampleMode::Argmax),
        &mut get_random_number_generator(0),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let mut sink = RecordingSeqSink::default();
      let derived = partition.stream_sequences(
        graph,
        SeqTrack::Nuc(1),
        true,
        report_unknown,
        Some(SeqStream::new(&mut sink, &names)),
      )?;
      let node_sequences: NodeMap<Seq> = sink.items.into_iter().map(|(key, _, seq)| (key, seq)).collect();
      let name = |key: GraphNodeKey| names[key].to_owned();
      let mut mutations = BTreeMap::new();
      for edge in graph.get_edges() {
        let subs = substitutions(&derived.edge_mutations[&edge.key()]);
        assert_rebuilds_child(
          &node_sequences[&edge.source()],
          &node_sequences[&edge.target()],
          &subs,
          &alphabet,
          report_unknown,
          &name(edge.target()),
        );
        mutations.insert(name(edge.target()), subs.iter().map(ToString::to_string).collect());
      }
      let sequences = node_sequences
        .iter()
        .map(|(key, seq)| (name(key), seq.as_str().to_owned()))
        .collect();
      Ok(Reconstructed { mutations, sequences })
    }

    fn assert_rebuilds_child(
      parent: &Seq,
      child: &Seq,
      subs: &[Sub],
      alphabet: &Alphabet,
      report_unknown: bool,
      child_name: &str,
    ) {
      let mut applied = parent.clone();
      for sub in subs {
        applied[sub.pos()] = sub.qry();
      }
      let child = child
        .iter()
        .zip(applied.iter())
        .map(|(&child, &applied)| {
          if alphabet.is_gap(applied) || (!report_unknown && child == alphabet.unknown()) {
            applied
          } else {
            child
          }
        })
        .collect::<Seq>();
      assert_eq!(
        child.as_str(),
        applied.as_str(),
        "parent sequence plus the mutations of the edge must give the sequence of {child_name}, with unknown tip states read as missing data unless reported, and gaps of the parent left to the indels"
      );
    }

    fn substitutions(mutations: &[Mutation]) -> Vec<Sub> {
      mutations
        .iter()
        .filter_map(|mutation| match &mutation.event {
          MutationEvent::Substitution(sub) => Some(sub.clone()),
          MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => None,
        })
        .collect()
    }

    pub(super) struct MaskedClade {
      pub masked_internal: BTreeMap<String, String>,
      pub sequences: BTreeMap<String, String>,
      pub mutations: BTreeMap<String, Vec<String>>,
    }

    pub(super) fn masked_clade_sequences(representation: Representation, impute: bool) -> Result<MaskedClade, Report> {
      let dataset = TestDataset::Flu20;
      let parse = nwk_read(dataset.tree().as_bytes())?;
      let clade = smallest_clade_with_two_internal_nodes(&parse.graph)?;
      let clade_leaves: BTreeSet<String> = clade.leaves.iter().map(|&key| parse.names[key].to_owned()).collect();
      let clade_internal: BTreeSet<String> = clade.internal.iter().map(|&key| parse.names[key].to_owned()).collect();
      let alphabet = Alphabet::default();
      let sequences = fasta_xz_records(dataset.alignment_xz(), &alphabet)?
        .into_iter()
        .map(|mut record| {
          if clade_leaves.contains(&record.seq_name) {
            record.seq[100..400].fill(alphabet.unknown());
          }
          AlignmentRecord::from(record)
        })
        .collect_vec();
      let Reconstructed { mutations, sequences } = reconstruct_records(
        dataset.tree(),
        sequences,
        ReconstructionMethod::Marginal(representation),
        &GtrModel::Infer,
        impute,
        true,
      )?;
      let masked_internal = sequences
        .iter()
        .filter(|(name, _)| clade_internal.contains(*name))
        .map(|(name, seq)| (name.clone(), seq.chars().skip(100).take(300).collect()))
        .collect();
      Ok(MaskedClade {
        masked_internal,
        sequences,
        mutations,
      })
    }

    struct Clade {
      leaves: Vec<GraphNodeKey>,
      internal: Vec<GraphNodeKey>,
    }

    fn smallest_clade_with_two_internal_nodes(graph: &Graph) -> Result<Clade, Report> {
      let schedule = GraphSchedule::new(graph)?;
      let root = graph.root_key()?;
      schedule
        .preorder()
        .iter()
        .filter(|&&key| key != root && !schedule.children(key).is_empty())
        .map(|&key| clade_below(&schedule, key))
        .filter(|clade| clade.internal.len() >= 2)
        .min_by_key(|clade| clade.leaves.len())
        .ok_or_else(|| make_report!("Tree has no clade with two internal nodes"))
    }

    fn clade_below(schedule: &GraphSchedule<'_>, key: GraphNodeKey) -> Clade {
      let mut leaves = vec![];
      let mut internal = vec![];
      let mut stack = vec![key];
      while let Some(key) = stack.pop() {
        let children = schedule.children(key);
        if children.is_empty() {
          leaves.push(key);
        } else {
          internal.push(key);
          stack.extend(children.iter().map(|&(child, _)| child));
        }
      }
      Clade { leaves, internal }
    }
  }
}
