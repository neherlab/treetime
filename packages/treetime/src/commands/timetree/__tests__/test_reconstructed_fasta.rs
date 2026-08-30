#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use treetime_primitives::AsciiChar;

  // Without a tip-state flag the reconstructed FASTA holds internal nodes only, and every record
  // carries a real name. Rerooting introduces a fresh root node, so this also guards that the root
  // reaches output named (`NODE_<n>`) rather than as an empty FASTA header.
  #[test]
  fn test_timetree_reconstructed_fasta_internal_only_and_named() -> Result<(), Report> {
    let records = helpers::run_reconstructed_fasta("test-tt-recon-internal", |_| {})?;
    let leaves = helpers::input_leaves()?;

    assert!(!records.is_empty(), "reconstruction must emit internal-node sequences");
    for record in &records {
      assert!(
        !record.seq_name.is_empty(),
        "every reconstructed record must have a non-empty name (the rerooted root included)"
      );
      assert!(
        !leaves.contains_key(&record.seq_name),
        "without --include-leaves the FASTA must not contain leaf {}",
        record.seq_name
      );
    }
    assert!(
      records.iter().any(|record| record.seq_name.starts_with("NODE_")),
      "internal nodes carry NODE_<n> names, matching the ancestral command and v0"
    );
    Ok(())
  }

  // `--include-leaves` adds every leaf, and without imputation each reconstructed leaf echoes its
  // observed input sequence exactly.
  #[test]
  fn test_timetree_reconstructed_fasta_include_leaves_echoes_observed() -> Result<(), Report> {
    let records = helpers::run_reconstructed_fasta("test-tt-recon-leaves", |args| args.include_leaves = true)?;
    let leaves = helpers::input_leaves()?;

    let emitted_leaf_names: BTreeSet<String> = records
      .iter()
      .filter(|record| leaves.contains_key(&record.seq_name))
      .map(|record| record.seq_name.clone())
      .collect();
    let expected_leaf_names: BTreeSet<String> = leaves.keys().cloned().collect();
    assert_eq!(
      expected_leaf_names, emitted_leaf_names,
      "--include-leaves must emit exactly the observed leaves"
    );

    for record in &records {
      if let Some(observed) = leaves.get(&record.seq_name) {
        assert_eq!(
          *observed, record.seq,
          "leaf {} must echo its observed sequence when imputation is off",
          record.seq_name
        );
      }
    }
    Ok(())
  }

  // `--reconstruct-tip-states` enables imputation: every ambiguous observed tip position resolves to
  // an inferred canonical state, while gaps stay gaps.
  #[test]
  fn test_timetree_reconstructed_fasta_imputes_ambiguous_tips() -> Result<(), Report> {
    let records = helpers::run_reconstructed_fasta("test-tt-recon-impute", |args| args.reconstruct_tip_states = true)?;
    let leaves = helpers::input_leaves()?;
    let alphabet = Alphabet::default();

    let mut resolved = 0_usize;
    for record in &records {
      let Some(observed) = leaves.get(&record.seq_name) else {
        continue;
      };
      let observed = observed.as_str().as_bytes();
      let imputed = record.seq.as_str().as_bytes();
      assert_eq!(observed.len(), imputed.len(), "imputation preserves sequence length");
      for (obs, imp) in observed.iter().zip(imputed) {
        let obs = AsciiChar::from_byte_unchecked(*obs);
        let imp = AsciiChar::from_byte_unchecked(*imp);
        if alphabet.is_gap(obs) {
          assert!(alphabet.is_gap(imp), "gaps are inferred structure and stay gaps");
        } else if alphabet.is_ambiguous(obs) || alphabet.is_unknown(obs) {
          assert!(
            !alphabet.is_ambiguous(imp) && !alphabet.is_unknown(imp),
            "imputed tip position must resolve to a canonical state, got {}",
            char::from(imp.inner())
          );
          resolved += 1;
        }
      }
    }
    assert!(
      resolved > 0,
      "the test dataset must contain ambiguous tip positions for imputation to resolve"
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::timetree::args::TreetimeTimetreeArgs;
    use crate::commands::timetree::run::run_timetree_estimation;
    use crate::progress::NoopProgress;
    use eyre::Report;
    use std::collections::BTreeMap;
    use std::path::PathBuf;
    use treetime_io::fasta::{FastaRecord, read_many_fasta};
    use treetime_primitives::Seq;

    pub fn project_root() -> PathBuf {
      PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .map(PathBuf::from)
        .expect("project has workspace root")
    }

    pub fn input_leaves() -> Result<BTreeMap<String, Seq>, Report> {
      let alignment = project_root().join("data/flu/h3n2/20/aln.fasta.xz");
      let records = read_many_fasta(&[alignment], &Alphabet::default())?;
      Ok(
        records
          .into_iter()
          .map(|record| (record.seq_name, record.seq))
          .collect(),
      )
    }

    pub fn run_reconstructed_fasta(
      subdir: &str,
      configure: impl FnOnce(&mut TreetimeTimetreeArgs),
    ) -> Result<Vec<FastaRecord>, Report> {
      let root = project_root();
      let outdir = root.join("tmp").join(subdir);
      std::fs::create_dir_all(&outdir)?;
      let fasta = outdir.join("ancestral_sequences.fasta");

      let mut args = TreetimeTimetreeArgs {
        alignment: AlignmentArgs {
          alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
        },
        tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
        metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
        max_iter: 2,
        output_reconstructed_nuc_fasta: Some(fasta.clone()),
        ..TreetimeTimetreeArgs::default()
      };
      configure(&mut args);

      run_timetree_estimation(&args, &NoopProgress)?;
      read_many_fasta(&[fasta], &Alphabet::default())
    }
  }
}
