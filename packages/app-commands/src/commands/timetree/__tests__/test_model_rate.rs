#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::project_root;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
  use crate::commands::shared::output_args::{OutputCoreArgs, TimetreeOutputSelection};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::run::run_timetree_estimation;
  use app_output::output_plan::OutputSelection;
  use app_output::table_output::table_read_file;
  use eyre::Report;
  use helpers::run_symmetric;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::Path;
  use treetime::cancel::NoopCancel;
  use treetime::clock::clock_model::ClockModel;
  use treetime::gtr::get_gtr::{ModelEntry, ModelOutput};
  use treetime::progress::NoopProgress;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_map_abs_diff_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( true )]
  #[case::sparse(false)]
  #[trace]
  fn test_timetree_one_alignment_with_inferred_model_has_unit_overall_rate(#[case] dense: bool) -> Result<(), Report> {
    let output = tempfile::tempdir()?;
    let dir = project_root().join("data/flu/h3n2/20");
    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![dir.join("aln.fasta.xz")],
      },
      tree: Some(dir.join("tree.nwk")),
      metadata: Some(dir.join("metadata.tsv")),
      dense: Some(dense),
      max_iter: 1,
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..OutputCoreArgs::default()
      },
      output_selection: vec![TimetreeOutputSelection::Gtr, TimetreeOutputSelection::ClockModel],
      ..TreetimeTimetreeArgsRaw::default()
    })?;

    run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let models: ModelOutput = json_read_file(output.path().join("timetree.gtr.json"))?;
    let clock_model: ClockModel = json_read_file(output.path().join("timetree.clock-model.json"))?;
    let [ModelEntry::Sequence(model)] = models.models.as_slice() else {
      panic!("one alignment gives one sequence model");
    };
    pretty_assert_ulps_eq!(
      1.0,
      model.model.mu(),
      max_ulps = 0,
      "v0 sets the overall rate of an inferred model to 1: packages/legacy/treetime/treetime/treeanc.py#L1626-L1631"
    );
    assert_eq!(Some(clock_model.clock_rate()), model.clock_rate);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( true )]
  #[case::sparse(false)]
  #[trace]
  fn test_timetree_inferred_model_on_nucleotide_symmetric_data_matches_jc69(#[case] dense: bool) -> Result<(), Report> {
    let inferred = run_symmetric(GtrModelNameCli::Infer, dense)?;
    let jc69 = run_symmetric(GtrModelNameCli::JC69, dense)?;

    pretty_assert_abs_diff_eq!(
      jc69.clock_rate,
      inferred.clock_rate,
      epsilon = 1e-13,
      "data symmetric under every permutation of ACGT infer JC69, whose branch lengths are substitutions per site (Jukes & Cantor 1969)"
    );
    pretty_assert_map_abs_diff_eq!(jc69.divergence, inferred.divergence, epsilon = 1e-13);
    Ok(())
  }

  mod helpers {
    use super::*;
    use std::iter::repeat_n;
    use treetime::clock::rtt::ClockRegressionResult;

    const SUBSTITUTIONS: [(&str, usize); 7] = [
      ("yxxxx", 1),
      ("xyxxx", 2),
      ("xxyxx", 2),
      ("xxxyx", 3),
      ("xxxxy", 1),
      ("yyxxx", 1),
      ("xxyyx", 1),
    ];
    const INVARIANT_COLUMNS_PER_STATE: usize = 100;
    const NUCLEOTIDES: [u8; 4] = *b"ACGT";

    pub(super) fn run_symmetric(model: GtrModelNameCli, dense: bool) -> Result<SymmetricRun, Report> {
      let dir = tempfile::tempdir()?;
      let alignment = dir.path().join("aln.fasta");
      let tree = dir.path().join("tree.nwk");
      let metadata = dir.path().join("metadata.tsv");
      fs::write(&alignment, symmetric_alignment())?;
      fs::write(
        &tree,
        "((A:0.0226,B:0.0451)AB:0.0226,(C:0.0451,D:0.0677)CD:0.0226,E:0.0226)R;\n",
      )?;
      fs::write(&metadata, "strain\tdate\nA\t2002\nB\t2003\nC\t2003\nD\t2004\nE\t2001\n")?;
      let output = dir.path().join("out");
      let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![alignment],
        },
        tree: Some(tree),
        metadata: Some(metadata),
        model_args: ModelArgs {
          model,
          ..ModelArgs::default()
        },
        dense: Some(dense),
        keep_root: true,
        output: OutputCoreArgs {
          output_all: Some(output.clone()),
          ..OutputCoreArgs::default()
        },
        output_selection: vec![TimetreeOutputSelection::ClockModel, TimetreeOutputSelection::ClockCsv],
        ..TreetimeTimetreeArgsRaw::default()
      })?;
      run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;
      read_run(&output)
    }

    pub(super) struct SymmetricRun {
      pub clock_rate: f64,
      pub divergence: BTreeMap<String, f64>,
    }

    fn symmetric_alignment() -> String {
      let pairs = NUCLEOTIDES.iter().permutations(2).collect_vec();
      let mut columns = NUCLEOTIDES
        .iter()
        .flat_map(|&state| repeat_n(vec![state; 5], INVARIANT_COLUMNS_PER_STATE))
        .collect_vec();
      for (pattern, count) in SUBSTITUTIONS {
        for pair in &pairs {
          let column = pattern
            .bytes()
            .map(|symbol| if symbol == b'x' { *pair[0] } else { *pair[1] })
            .collect_vec();
          columns.extend(repeat_n(column, count));
        }
      }
      ["A", "B", "C", "D", "E"]
        .iter()
        .enumerate()
        .map(|(row, name)| {
          let sequence: String = columns.iter().map(|column| char::from(column[row])).collect();
          format!(">{name}\n{sequence}\n")
        })
        .join("")
    }

    fn read_run(output: &Path) -> Result<SymmetricRun, Report> {
      let clock_model: ClockModel = json_read_file(output.join("timetree.clock-model.json"))?;
      let rows: Vec<ClockRegressionResult> =
        table_read_file(OutputSelection::ClockCsv, &output.join("timetree.clock.csv"))?;
      let divergence = rows
        .into_iter()
        .filter_map(|row| row.name.map(|name| (name, row.div)))
        .collect();
      Ok(SymmetricRun {
        clock_rate: clock_model.clock_rate(),
        divergence,
      })
    }
  }
}
