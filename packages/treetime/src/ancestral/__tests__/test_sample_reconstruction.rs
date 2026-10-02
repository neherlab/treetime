#[cfg(test)]
mod tests {
  use crate::ancestral::sample::SampleMode;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_sample_reconstruction_argmax_ignores_seed() -> Result<(), Report> {
    let run_a = helpers::reconstruct(SampleMode::Argmax, 1)?;
    let run_b = helpers::reconstruct(SampleMode::Argmax, 2)?;
    assert_eq!(run_a, run_b);
    Ok(())
  }

  #[test]
  fn test_sample_reconstruction_all_seeded_reproducible() -> Result<(), Report> {
    let run_a = helpers::reconstruct(SampleMode::All, 12345)?;
    let run_b = helpers::reconstruct(SampleMode::All, 12345)?;
    assert_eq!(run_a, run_b);
    Ok(())
  }

  #[test]
  fn test_sample_reconstruction_root_seeded_reproducible() -> Result<(), Report> {
    let run_a = helpers::reconstruct(SampleMode::Root, 777)?;
    let run_b = helpers::reconstruct(SampleMode::Root, 777)?;
    assert_eq!(run_a, run_b);
    Ok(())
  }

  #[test]
  fn test_sample_reconstruction_root_only_leaves_nonroot_unchanged() -> Result<(), Report> {
    let argmax = helpers::reconstruct(SampleMode::Argmax, 0)?;
    let root_sampled = helpers::reconstruct(SampleMode::Root, 777)?;

    assert_eq!(
      argmax.keys().collect::<Vec<_>>(),
      root_sampled.keys().collect::<Vec<_>>()
    );

    for (name, seq) in &argmax {
      if name == helpers::ROOT_NAME {
        continue;
      }
      assert_eq!(
        seq, &root_sampled[name],
        "non-root node {name} must match argmax under root sampling"
      );
    }
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::fitch::create_fitch_partition;
    use crate::ancestral::marginal::branch_lengths_or_zero;
    use crate::ancestral::pipeline::SparseReconstruction;
    use crate::ancestral::sample::SampleMode;
    use crate::ancestral::tip_states::TipStates;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::seq::alignment::node_seq_inputs;
    use crate::test_utils::emitted_sequences_by_name;
    use eyre::Report;
    use indoc::indoc;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use std::collections::BTreeMap;
    use treetime_graph::graph::Graph;
    use treetime_io::fasta::read_many_fasta_str;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::AlignmentRecord;

    pub(super) const ROOT_NAME: &str = "root";

    const TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    pub(super) fn reconstruct(mode: SampleMode, seed: u64) -> Result<BTreeMap<String, String>, Report> {
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
        &Alphabet::default(),
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

      let nwk_parsed = nwk_read_str(TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let graph: Graph = graph;
      let fitch = create_fitch_partition(&graph, 0, Alphabet::default(), &node_seq_inputs(&graph, &names, aln))?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let recon = SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);

      let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

      let reconstruction =
        recon.reconstruct_sequences(&graph, TipStates::default(), mode, &mut StdRng::seed_from_u64(seed))?;
      Ok(
        emitted_sequences_by_name(&names, &reconstruction, |key| recon.node_sequence(&graph, false, key))?
          .into_iter()
          .map(|(name, seq)| (name, seq.to_string()))
          .collect(),
      )
    }
  }
}
