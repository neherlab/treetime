#[cfg(test)]
mod tests {
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output_args::{OutputCoreArgs, TimetreeOutputSelection};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::run::run_timetree_estimation;
  use approx::{assert_relative_eq, assert_ulps_eq};
  use eyre::Report;
  use helpers::{fixed_rate_intercept, ordinary_least_squares, run_zika, zika_dir};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::fs;
  use std::path::{Path, PathBuf};
  use treetime::cancel::NoopCancel;
  use treetime::clock::clock_model::ClockModel;
  use treetime::clock::rtt::{ClockDateSource, ClockRegressionResult};
  use treetime::o;
  use treetime::progress::NoopProgress;
  use treetime_io::csv::{csv_read_file, default_metadata_delimiters, default_name_candidates};
  use treetime_io::dates_csv::read_dates;
  use treetime_io::nwk::nwk_read_file;
  use treetime_utils::io::json::json_read_file;

  const UNDATED_TIP: &str = "GZ01|KU820898|2016-02-14|china";

  #[rustfmt::skip]
  #[rstest]
  #[case::refined(         (2, false, None))]
  #[case::rerooted_only(   (0, false, None))]
  #[case::input_root_only( (0, true,  None))]
  #[case::fixed_rate(      (2, false, Some(8e-4)))]
  #[trace]
  fn test_timetree_clock_csv_refits_the_clock_model(
    #[case] (max_iter, keep_root, clock_rate): (usize, bool, Option<f64>),
  ) -> Result<(), Report> {
    let output = tempfile::tempdir()?;
    let (rows, model) = run_zika(output.path(), &zika_dir().join("metadata.tsv"), |raw| {
      raw.max_iter = max_iter;
      raw.keep_root = keep_root;
      raw.clock_rate = clock_rate;
    })?;

    let points = rows
      .iter()
      .filter(|row| !row.is_outlier)
      .filter_map(|row| row.date.map(|date| (date, row.div)))
      .collect_vec();
    let (rate, intercept) = match clock_rate {
      Some(rate) => (rate, fixed_rate_intercept(&points, rate)),
      None => ordinary_least_squares(&points),
    };

    assert_relative_eq!(model.clock_rate(), rate, max_relative = 1e-8);
    assert_relative_eq!(model.intercept(), intercept, max_relative = 1e-8);
    Ok(())
  }

  #[test]
  fn test_timetree_clock_csv_lists_every_tip_with_its_metadata_date() -> Result<(), Report> {
    let output = tempfile::tempdir()?;
    let metadata = zika_dir().join("metadata.tsv");
    let (rows, _) = run_zika(output.path(), &metadata, |_| {})?;

    let tree = nwk_read_file(zika_dir().join("tree.nwk"))?;
    let names = tree.names();
    let tips: BTreeSet<String> = tree
      .graph
      .get_leaves()
      .filter_map(|leaf| names[&leaf.key()].clone())
      .collect();
    let dates = read_dates(
      &metadata,
      &default_metadata_delimiters(),
      &default_name_candidates(),
      &None,
      &Some(o!("date")),
    )?;

    let listed: BTreeSet<String> = rows.iter().filter_map(|row| row.name.clone()).collect();
    assert_eq!(20, rows.len());
    assert_eq!(tips, listed);
    assert_eq!(tips, dates.keys().cloned().collect::<BTreeSet<_>>());
    for row in &rows {
      let name = row.name.as_deref().expect("every tip has a name");
      let expected = dates[name].as_ref().expect("zika/20 dates every sample").mean();
      assert_eq!(Some(ClockDateSource::Input), row.date_source, "{name}");
      assert_ulps_eq!(expected, row.date.expect("an input date"), max_ulps = 4);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::refined(       (2, Some(ClockDateSource::Inferred), true))]
  #[case::rerooted_only( (0, Some(ClockDateSource::Missing),  false))]
  #[trace]
  fn test_timetree_clock_csv_marks_the_date_of_an_undated_tip(
    #[case] (max_iter, source, dated): (usize, Option<ClockDateSource>, bool),
  ) -> Result<(), Report> {
    let output = tempfile::tempdir()?;
    let metadata = output.path().join("metadata.tsv");
    let text = fs::read_to_string(zika_dir().join("metadata.tsv"))?;
    let undated = text
      .lines()
      .map(|line| {
        if line.starts_with(UNDATED_TIP) {
          format!("{UNDATED_TIP}\t\tchina")
        } else {
          line.to_owned()
        }
      })
      .join("\n");
    fs::write(&metadata, undated)?;

    let (rows, _) = run_zika(&output.path().join("out"), &metadata, |raw| raw.max_iter = max_iter)?;

    let row = rows
      .iter()
      .find(|row| row.name.as_deref() == Some(UNDATED_TIP))
      .expect("the undated tip is listed");
    assert_eq!((source, dated), (row.date_source, row.date.is_some()));
    assert!(
      rows
        .iter()
        .filter(|row| row.name.as_deref() != Some(UNDATED_TIP))
        .all(|row| row.date_source == Some(ClockDateSource::Input))
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn run_zika(
      output: &Path,
      metadata: &Path,
      configure: impl FnOnce(&mut TreetimeTimetreeArgsRaw),
    ) -> Result<(Vec<ClockRegressionResult>, ClockModel), Report> {
      let mut raw = TreetimeTimetreeArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![zika_dir().join("aln.fasta.xz")],
        },
        tree: Some(zika_dir().join("tree.nwk")),
        metadata: Some(metadata.to_path_buf()),
        output: OutputCoreArgs {
          output_all: Some(output.to_path_buf()),
          ..OutputCoreArgs::default()
        },
        output_selection: vec![TimetreeOutputSelection::ClockModel, TimetreeOutputSelection::ClockCsv],
        seed: Some(7),
        ..TreetimeTimetreeArgsRaw::default()
      };
      configure(&mut raw);
      let args = TreetimeTimetreeArgs::try_from(raw)?;
      run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;
      let rows = csv_read_file(output.join("timetree.clock.csv"), b',')?;
      let model = json_read_file(output.join("timetree.clock-model.json"))?;
      Ok((rows, model))
    }

    pub(super) fn ordinary_least_squares(points: &[(f64, f64)]) -> (f64, f64) {
      let (mean_t, mean_d) = means(points);
      let sxy: f64 = points.iter().map(|(t, d)| (t - mean_t) * (d - mean_d)).sum();
      let sxx: f64 = points.iter().map(|(t, _)| (t - mean_t).powi(2)).sum();
      let rate = sxy / sxx;
      (rate, mean_d - rate * mean_t)
    }

    pub(super) fn fixed_rate_intercept(points: &[(f64, f64)], rate: f64) -> f64 {
      let (mean_t, mean_d) = means(points);
      mean_d - rate * mean_t
    }

    #[allow(clippy::as_conversions, reason = "point count is far below 2^52")]
    pub(super) fn means(points: &[(f64, f64)]) -> (f64, f64) {
      let n = points.len() as f64;
      (
        points.iter().map(|(t, _)| t).sum::<f64>() / n,
        points.iter().map(|(_, d)| d).sum::<f64>() / n,
      )
    }

    pub(super) fn zika_dir() -> PathBuf {
      PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20")
    }
  }
}
