#[cfg(test)]
mod tests {
  use crate::ancestral::sample::SampleMode;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  /// Argmax reconstruction ignores the RNG entirely: different seeds yield identical sequences.
  #[test]
  fn test_sample_reconstruction_argmax_ignores_seed() -> Result<(), Report> {
    let run_a = helpers::reconstruct(SampleMode::Argmax, 1)?;
    let run_b = helpers::reconstruct(SampleMode::Argmax, 2)?;
    assert_eq!(run_a, run_b);
    Ok(())
  }

  /// Sampling at every node is reproducible: the same seed reproduces the same draws and sequences.
  #[test]
  fn test_sample_reconstruction_all_seeded_reproducible() -> Result<(), Report> {
    let run_a = helpers::reconstruct(SampleMode::All, 12345)?;
    let run_b = helpers::reconstruct(SampleMode::All, 12345)?;
    assert_eq!(run_a, run_b);
    Ok(())
  }

  /// Root-only sampling is reproducible under a fixed seed.
  #[test]
  fn test_sample_reconstruction_root_seeded_reproducible() -> Result<(), Report> {
    let run_a = helpers::reconstruct(SampleMode::Root, 777)?;
    let run_b = helpers::reconstruct(SampleMode::Root, 777)?;
    assert_eq!(run_a, run_b);
    Ok(())
  }

  /// Root-only sampling perturbs at most the root: every non-root internal node is identical to the
  /// deterministic argmax reconstruction.
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
    use crate::ancestral::marginal::{ancestral_reconstruction_marginal, update_marginal};
    use crate::ancestral::sample::SampleMode;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::payload::ancestral::GraphAncestral;
    use eyre::Report;
    use indoc::indoc;
    use parking_lot::RwLock;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use std::collections::BTreeMap;
    use std::sync::Arc;
    use treetime_io::fasta::read_many_fasta_str;
    use treetime_io::nwk::nwk_read_str;

    pub const ROOT_NAME: &str = "root";

    const TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    /// Build a fresh sparse marginal partition, run the marginal passes, and reconstruct internal
    /// node sequences with the given sampling mode and seed. Returns a name -> sequence map.
    ///
    /// A fresh partition is built on every call so that the seeded RNG is the only source of
    /// variation between runs, making reproducibility assertions meaningful.
    pub fn reconstruct(mode: SampleMode, seed: u64) -> Result<BTreeMap<String, String>, Report> {
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
        &Alphabet::default(),
      )?;

      let graph: GraphAncestral = nwk_read_str(TREE)?;
      let fitch = create_fitch_partition(&graph, 0, Alphabet::default(), &aln)?;
      let partitions = [Arc::new(RwLock::new(
        fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
      ))];

      update_marginal(&graph, &partitions)?;

      let mut rng = StdRng::seed_from_u64(seed);
      let mut out = BTreeMap::new();
      ancestral_reconstruction_marginal(&graph, false, &partitions, mode, &mut rng, |node, seq| {
        out.insert(node.name.clone().unwrap_or_default(), seq.to_string());
        Ok(())
      })?;
      Ok(out)
    }
  }
}
