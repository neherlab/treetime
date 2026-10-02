#[cfg(test)]
mod tests {
  use crate::partition::marginal::sample::SampleMode;
  use eyre::{OptionExt, Report};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;

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

  #[test]
  fn test_sample_reconstruction_sparse_constant_sites_always_sample_observed_state() -> Result<(), Report> {
    for seed in 0..64 {
      let root = helpers::reconstruct_two_leaves(SampleMode::Root, seed)?;
      assert_eq!(
        "ACGT",
        root.chars().take(4).collect::<String>(),
        "constant columns are certain at seed {seed}"
      );
    }
    Ok(())
  }

  #[test]
  fn test_sample_reconstruction_two_state_site_samples_both_states() -> Result<(), Report> {
    let sampled: BTreeSet<char> = (0..64)
      .map(|seed| {
        let root = helpers::reconstruct_two_leaves(SampleMode::Root, seed)?;
        root.chars().nth(4).ok_or_eyre("root sequence must have five sites")
      })
      .collect::<Result<_, Report>>()?;
    assert!(
      sampled.contains(&'A') && sampled.contains(&'C'),
      "the variable column has equal posterior on A and C, sampled states: {sampled:?}"
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::reconstruction::{MarginalReconstruction, SparseReconstruction};
    use crate::partition::marginal::sample::SampleMode;
    use crate::partition::marginal::sequences::TipStates;
    use crate::seq::alignment::node_seq_inputs;
    use crate::test_utils::emitted_sequences_by_name;
    use eyre::{OptionExt, Report};
    use indoc::indoc;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use std::collections::BTreeMap;
    use treetime_graph::graph::Graph;
    use treetime_io::fasta::read_many_fasta_str;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::AlignmentRecord;

    pub(super) const ROOT_NAME: &str = "root";

    pub(super) fn reconstruct(mode: SampleMode, seed: u64) -> Result<BTreeMap<String, String>, Report> {
      reconstruct_named(
        "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;",
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
        mode,
        seed,
      )
    }

    pub(super) fn reconstruct_two_leaves(mode: SampleMode, seed: u64) -> Result<String, Report> {
      let mut sequences = reconstruct_named(
        "(A:0.1,B:0.1)root;",
        indoc! {r#"
        >A
        ACGTA
        >B
        ACGTC
      "#},
        mode,
        seed,
      )?;
      sequences.remove(ROOT_NAME).ok_or_eyre("root must be emitted")
    }

    fn reconstruct_named(
      tree: &str,
      fasta: &str,
      mode: SampleMode,
      seed: u64,
    ) -> Result<BTreeMap<String, String>, Report> {
      let aln: Vec<AlignmentRecord> = read_many_fasta_str(fasta, &Alphabet::default())?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();

      let nwk_parsed = nwk_read_str(tree)?;
      let names = nwk_parsed.names();
      let graph: Graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let fitch = create_fitch_partition(&graph, 0, Alphabet::default(), &node_seq_inputs(&graph, &names, aln))?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let recon = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
        partition,
        jc69(JC69Params::default())?,
        node_states,
      ));

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
