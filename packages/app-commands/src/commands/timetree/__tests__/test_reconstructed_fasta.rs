#[cfg(test)]
mod tests {
  use crate::commands::shared::argument_error::ArgumentError;
  use crate::commands::timetree::branch_length_mode::{
    BranchLengthModeCli, ResolvedBranchLengthMode, settle_branch_length_outputs,
  };
  use app_output::output_plan::{OutputSelection, Requested};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use std::path::PathBuf;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::progress::NoopProgress;
  use treetime_primitives::AsciiChar;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_timetree_reconstructed_fasta_requires_ancestral_reconstruction() {
    let result = helpers::run_reconstructed_fasta("test-tt-recon-input-mode", |args| {
      args.branch_length_mode = BranchLengthModeCli::Input;
    });
    assert_error!(
      result,
      "--output-sequences was requested, but --branch-length-mode input keeps the input branch lengths and reconstructs \
       no sequences"
    );
  }

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

  #[test]
  fn test_timetree_reconstructed_fasta_include_leaves_echoes_observed() -> Result<(), Report> {
    let records = helpers::run_reconstructed_fasta("test-tt-recon-leaves", |args| {
      args.tip_sequences.include_leaves = true;
    })?;
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

  #[test]
  fn test_timetree_reconstructed_fasta_imputes_ambiguous_tips() -> Result<(), Report> {
    let records = helpers::run_reconstructed_fasta("test-tt-recon-impute", |args| {
      args.tip_sequences.include_leaves = true;
      args.impute_missing_data.impute_missing_data = true;
    })?;
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
        } else if alphabet.is_ambiguous(obs) || obs == alphabet.unknown() {
          assert!(
            !alphabet.is_ambiguous(imp) && imp != alphabet.unknown(),
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

  #[rustfmt::skip]
  #[rstest]
  #[case::marginal_named((Some(0.05), Requested::Named), btreemap! { 1 => PathBuf::from("out.fasta") })]
  #[case::marginal_all(  (Some(0.05), Requested::All),   btreemap! { 1 => PathBuf::from("out.fasta") })]
  #[case::no_length_all( (None,       Requested::All),   btreemap! { 1 => PathBuf::from("out.fasta") })]
  #[case::input_all(     (Some(0.25), Requested::All),   btreemap! {})]
  #[trace]
  fn test_timetree_conditional_fasta_follows_the_resolved_mode(
    #[case] (longest_branch, requested): (Option<f64>, Requested),
    #[case] expected: BTreeMap<usize, PathBuf>,
  ) -> Result<(), Report> {
    let mut resolved = helpers::conditional_fasta_outputs(requested);
    let mode = ResolvedBranchLengthMode::new(BranchLengthModeCli::Auto, true, longest_branch);

    settle_branch_length_outputs(&mut resolved, &mode, &NoopProgress)?;

    assert_eq!(
      (expected, true),
      (resolved.alignment_paths(OutputSelection::Sequences), resolved.conditional_selections().is_empty())
    );
    Ok(())
  }

  #[test]
  fn test_timetree_conditional_fasta_named_fails_when_auto_resolves_to_input() {
    let mut resolved = helpers::conditional_fasta_outputs(Requested::Named);
    let mode = ResolvedBranchLengthMode::new(BranchLengthModeCli::Auto, true, Some(0.25));

    let result = settle_branch_length_outputs(&mut resolved, &mode, &NoopProgress);

    let Err(ArgumentError { message, settings }) = result.map_err(|report| report.downcast::<ArgumentError>().unwrap())
    else {
      panic!("expected an argument error");
    };
    assert_eq!(
      (
        "--output-sequences requires re-estimated branch lengths (--branch-length-mode marginal), but \
         --branch-length-mode auto chose input: the longest input branch is 2.500e-1 substitutions per site",
        vec![o!("output_sequences"), o!("branch_length_mode")]
      ),
      (message.as_str(), settings)
    );
  }

  mod helpers {
    use crate::__tests__::test_support::tests::project_root;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::metadata::MetadataArgs;
    use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
    use crate::commands::shared::seed::SeedArgs;
    use crate::commands::shared::tree_input::TreeArgs;
    use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
    use crate::commands::timetree::run::run_timetree_estimation;
    use app_output::output_conditions::OutputCondition;
    use app_output::output_plan::{OutputSelection, PlannedFile, Requested, ResolvedOutputs};
    use eyre::{Report, WrapErr};
    use maplit::btreemap;
    use std::collections::BTreeMap;
    use std::num::NonZeroUsize;
    use std::path::PathBuf;
    use treetime::alphabet::alphabet::Alphabet;
    use treetime::cancel::NoopCancel;
    use treetime::progress::NoopProgress;
    use treetime_io::fasta::{FastaRecord, fasta_read_file};
    use treetime_primitives::Seq;

    pub(super) fn input_leaves() -> Result<BTreeMap<String, Seq>, Report> {
      let alignment = project_root().join("data/flu/h3n2/20/aln.fasta.xz");
      let records = fasta_read_file(alignment, &Alphabet::default())?;
      Ok(
        records
          .into_iter()
          .map(|record| (record.seq_name, record.seq))
          .collect(),
      )
    }

    pub(super) fn run_reconstructed_fasta(
      subdir: &str,
      configure: impl FnOnce(&mut TreetimeTimetreeArgsRaw),
    ) -> Result<Vec<FastaRecord>, Report> {
      let root = project_root();
      let outdir = root.join("tmp").join(subdir);
      std::fs::create_dir_all(&outdir).wrap_err("When creating the output directory")?;
      let fasta = outdir.join("ancestral_sequences.fasta");

      let mut raw = TreetimeTimetreeArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
        },
        tree: TreeArgs {
          tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
        },
        metadata: MetadataArgs {
          metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
        },
        max_iter: 2,
        reconstructed_fasta_output: ReconstructedFastaOutputArgs {
          output_sequences: Some(fasta.clone()),
        },
        seed_args: SeedArgs { seed: Some(7) },
        ..TreetimeTimetreeArgsRaw::default()
      };
      configure(&mut raw);
      let args = TreetimeTimetreeArgs::try_from(raw)?;

      run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;
      fasta_read_file(fasta, &Alphabet::default())
    }

    pub(super) fn conditional_fasta_outputs(requested: Requested) -> ResolvedOutputs {
      ResolvedOutputs {
        alignment_count: NonZeroUsize::MIN,
        tree_outputs: btreemap! {},
        non_tree_outputs: btreemap! {},
        tree_outputs_by_alignment: btreemap! {},
        non_tree_outputs_by_alignment: btreemap! {
          OutputSelection::Sequences => vec![PlannedFile::new(
            PathBuf::from("out.fasta"),
            requested,
            Some(OutputCondition::BranchLengthsReestimated),
          )],
        },
      }
    }
  }
}
