#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::json_float::JsonFloat;
  use crate::results::auspice::run_auspice;
  use crate::results::clades::{CladeRequest, clade_in_runs};
  use crate::results::clock::{ClockLine, RootToTip};
  use crate::results::compare::compare_runs;
  use crate::results::run_results::{CommandResults, RunResults, run_results};
  use crate::runs::headline::RunHeadline;
  use crate::runs::manager::RunManager;
  use crate::runs::record::CreateRunRequest;
  use crate::runs::setting_differences::SettingDifference;
  use eyre::Report;
  use helpers::{auspice_path, clock_model, dataset, finished_run, timetree_config, zika};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::fs;
  use tempfile::tempdir;
  use treetime::clock::rtt::ClockDateSource;
  use treetime::o;
  use treetime_utils::assert_error;

  #[test]
  fn test_run_results_of_a_timetree_run_hold_estimates_and_regression() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, AppCommand::Timetree, timetree_config(None));

    let results = run_results(&runs, &id)?;

    let model = clock_model(&runs, &id)?;
    let RunResults {
      tree,
      results: CommandResults::Timetree(timetree),
      problems,
      ..
    } = results
    else {
      panic!("a timetree run has timetree results");
    };
    let estimates = timetree.estimates.expect("the run wrote an Auspice tree");
    let RootToTip { points, line } = timetree.root_to_tip.expect("app runs write the clock regression table");
    assert_eq!(
      (
        Vec::<String>::new(),
        20,
        20,
        20,
        Some(ClockLine {
          rate: model.clock_rate(),
          intercept: model.intercept()
        }),
        Some(model.clock_rate()),
      ),
      (
        problems.into_iter().map(|problem| problem.message).collect(),
        tree.expect("the run wrote an Auspice tree").tips().count(),
        estimates.samples,
        points
          .iter()
          .filter(|point| point.date_source == Some(ClockDateSource::Input))
          .count(),
        line,
        estimates.clock_rate,
      )
    );
    assert_eq!(timetree.iterations.len(), estimates.iterations);
    assert!(estimates.log_likelihood.is_some());
    Ok(())
  }

  #[test]
  fn test_run_results_headline_repeats_the_estimates_of_the_results() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, AppCommand::Timetree, timetree_config(None));

    let CommandResults::Timetree(timetree) = run_results(&runs, &id)?.results else {
      panic!("a timetree run has timetree results");
    };
    let estimates = timetree.estimates.expect("the run wrote an Auspice tree");
    assert_eq!(
      RunHeadline {
        root_date: estimates.root_date,
        clock_rate: estimates.clock_rate.map(JsonFloat),
        r_squared: estimates.r_squared.map(JsonFloat),
      },
      runs.get(&id)?.headline
    );
    Ok(())
  }

  #[test]
  fn test_run_results_of_a_clock_run_count_dated_samples_and_outliers() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(
      &runs,
      AppCommand::Clock,
      json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv") }),
    );

    let CommandResults::Clock(clock) = run_results(&runs, &id)?.results else {
      panic!("a clock run has clock results");
    };

    let points = clock
      .root_to_tip
      .expect("the clock command writes its regression table")
      .points;
    assert_eq!(
      (20, points.iter().filter(|point| point.outlier).count(), 20),
      (clock.estimates.dated_samples, clock.estimates.outliers, points.len())
    );
    Ok(())
  }

  #[test]
  fn test_run_results_of_a_mugration_run_count_the_states_of_the_samples() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(
      &runs,
      AppCommand::Mugration,
      json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv"), "attribute": "country" }),
    );

    let CommandResults::Mugration(mugration) = run_results(&runs, &id)?.results else {
      panic!("a mugration run has mugration results");
    };

    assert_eq!(
      (
        "country",
        10,
        mugration.state_changes.iter().map(|change| change.branches).sum()
      ),
      (
        mugration.attribute.as_str(),
        mugration.states,
        mugration.changed_branches
      )
    );
    assert!(mugration.root.is_some());
    Ok(())
  }

  #[test]
  fn test_run_results_of_an_ancestral_run_sum_the_mutations_of_the_branches() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(
      &runs,
      AppCommand::Ancestral,
      json!({ "tree": zika("tree.nwk"), "alignment": [zika("aln.fasta.xz")] }),
    );

    let results = run_results(&runs, &id)?;
    let CommandResults::Ancestral(ancestral) = &results.results else {
      panic!("an ancestral run has ancestral results");
    };

    let on_tree: usize = results
      .tree
      .as_ref()
      .map_or(0, |tree| tree.nodes.iter().map(|node| node.mutations.len()).sum());
    assert_eq!(
      (on_tree, true, true),
      (
        ancestral.mutations,
        ancestral.mutations > 0,
        ancestral.recurrent_sites.iter().all(|site| site.branches > 1)
      )
    );
    Ok(())
  }

  #[test]
  fn test_run_results_refuse_a_run_that_did_not_finish() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let created = runs.create(CreateRunRequest {
      command: AppCommand::Clock,
      config: json!({ "tree": zika("tree.nwk") }),
      title: None,
      defer_start: true,
    })?;

    assert_error!(
      run_results(&runs, &created.id),
      format!("run `{}` has no results: it is created", created.id.as_str())
    );
    Ok(())
  }

  #[test]
  fn test_run_auspice_colors_a_mugration_trait_with_the_app_palette() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(
      &runs,
      AppCommand::Mugration,
      json!({
        "tree": dataset("dengue/20", "tree.nwk"),
        "metadata": dataset("dengue/20", "metadata.tsv"),
        "metadata_id_columns": ["genbank_accession"],
        "attribute": "region",
      }),
    );

    let document = serde_json::to_value(run_auspice(&runs, &id)?)?;

    let tree = run_results(&runs, &id)?.tree.expect("the run wrote an Auspice tree");
    let region = tree
      .colorings
      .iter()
      .find(|coloring| coloring.key == "region")
      .expect("the trait is a coloring");
    let scale = document["meta"]["colorings"]
      .as_array()
      .expect("colorings")
      .iter()
      .find(|coloring| coloring["key"] == "region")
      .map(|coloring| coloring["scale"].clone());
    assert_eq!(
      (
        vec!["#332288", "#88ccee", "#44aa99", "#117733"],
        Some(json!(
          region
            .scale
            .iter()
            .map(|color| [&color.state, &color.color])
            .collect::<Vec<_>>()
        )),
      ),
      (
        region
          .scale
          .iter()
          .map(|color| color.color.as_str())
          .collect::<Vec<_>>(),
        scale,
      )
    );
    assert_eq!(
      region.states,
      region.scale.iter().map(|color| color.state.clone()).collect::<Vec<_>>()
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::mugration(AppCommand::Mugration, json!({ "tree": dataset("dengue/20", "tree.nwk"), "metadata": dataset("dengue/20", "metadata.tsv"), "metadata_id_columns": ["genbank_accession"], "attribute": "region" }))]
  #[case::timetree( AppCommand::Timetree,  timetree_config(None))]
  #[case::ancestral(AppCommand::Ancestral, json!({ "tree": zika("tree.nwk"), "alignment": [zika("aln.fasta.xz")] }))]
  #[case::clock(    AppCommand::Clock,     json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv") }))]
  #[trace]
  fn test_run_auspice_changes_only_color_scales_and_leaves_the_written_file_unchanged(
    #[case] command: AppCommand,
    #[case] config: Value,
  ) -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, command, config);
    let path = auspice_path(&runs, &id);
    let written = fs::read(&path)?;

    let document = serde_json::to_value(run_auspice(&runs, &id)?)?;

    let colorings = run_results(&runs, &id)?.tree.expect("the run wrote an Auspice tree").colorings;
    let mut expected: Value = serde_json::from_slice(&written)?;
    for (coloring, result) in expected["meta"]["colorings"].as_array_mut().expect("colorings").iter_mut().zip(&colorings) {
      if !result.scale.is_empty() {
        coloring["scale"] = json!(result.scale.iter().map(|color| [&color.state, &color.color]).collect::<Vec<_>>());
      }
    }
    assert_eq!((written, expected), (fs::read(&path)?, document));
    Ok(())
  }

  #[test]
  fn test_run_auspice_refuses_a_run_that_did_not_finish() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let created = runs.create(CreateRunRequest {
      command: AppCommand::Clock,
      config: json!({ "tree": zika("tree.nwk") }),
      title: None,
      defer_start: true,
    })?;

    assert_error!(
      run_auspice(&runs, &created.id),
      format!("run `{}` has no results: it is created", created.id.as_str())
    );
    Ok(())
  }

  #[test]
  fn test_compare_a_run_with_itself_shifts_nothing() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, AppCommand::Timetree, timetree_config(None));

    let comparison = compare_runs(&runs, &id, &id)?;

    let settings = comparison.settings.expect("both runs execute the same command");
    let estimates = comparison.estimates.expect("both runs are time trees");
    let ancestors = comparison.ancestors.expect("both runs are time trees");
    assert_eq!(
      (true, true),
      (settings.differences.is_empty(), settings.same_config_hash)
    );
    assert_eq!(
      (Some(0.0), Some(0.0), 0, ancestors.ancestors, true),
      (
        estimates.root_shift_days,
        estimates.clock_rate_change_percent,
        estimates.excluded_samples_change,
        ancestors.shifts.len(),
        ancestors.shifts.iter().all(|shift| shift.shift_days == 0.0)
      )
    );
    Ok(())
  }

  #[test]
  fn test_compare_runs_lists_the_settings_that_differ() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let first = finished_run(&runs, AppCommand::Timetree, timetree_config(None));
    let second = finished_run(&runs, AppCommand::Timetree, timetree_config(Some(1e-3)));

    let settings = compare_runs(&runs, &first, &second)?
      .settings
      .expect("both runs execute the same command");

    assert_eq!(
      (
        vec![SettingDifference::Setting {
          key: o!("clock_rate"),
          first: Value::Null,
          second: json!(1e-3),
        }],
        false,
        true,
      ),
      (settings.differences, settings.same_config_hash, settings.compared > 1)
    );
    Ok(())
  }

  #[test]
  fn test_compare_runs_of_other_commands_has_no_estimates() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let timetree = finished_run(&runs, AppCommand::Timetree, timetree_config(None));
    let clock = finished_run(
      &runs,
      AppCommand::Clock,
      json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv") }),
    );

    let comparison = compare_runs(&runs, &timetree, &clock)?;

    assert_eq!(
      (None, None, None),
      (comparison.settings, comparison.estimates, comparison.ancestors)
    );
    Ok(())
  }

  #[test]
  fn test_clade_in_runs_finds_the_root_clade_in_other_time_trees() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let first = finished_run(&runs, AppCommand::Timetree, timetree_config(None));
    let second = finished_run(&runs, AppCommand::Timetree, timetree_config(Some(1e-3)));
    finished_run(
      &runs,
      AppCommand::Clock,
      json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv") }),
    );
    let first_tree = run_results(&runs, &first)?.tree.expect("a tree");
    let second_tree = run_results(&runs, &second)?.tree.expect("a tree");

    let found = clade_in_runs(
      &runs,
      &CladeRequest {
        run: first,
        node: first_tree.root().name.clone(),
      },
    )?;

    assert_eq!(
      (
        1,
        vec![(second, second_tree.root().name.clone(), second_tree.root().date.clone())]
      ),
      (
        found.searched_runs,
        found
          .matches
          .into_iter()
          .map(|found| (found.run, found.node, found.date))
          .collect::<Vec<_>>()
      )
    );
    Ok(())
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::job::JobId;
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, RunStatus};
    use app_output::output_plan::OutputSelection;
    use eyre::Report;
    use serde_json::{Value, json};
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use treetime::clock::clock_model::ClockModel;
    use treetime_utils::io::json::json_read_file;

    pub(super) fn timetree_config(clock_rate: Option<f64>) -> Value {
      json!({
        "tree": zika("tree.nwk"),
        "metadata": zika("metadata.tsv"),
        "alignment": [zika("aln.fasta.xz")],
        "max_iter": 2,
        "seed": 7,
        "clock_rate": clock_rate,
      })
    }

    pub(super) fn clock_model(runs: &Arc<RunManager>, id: &JobId) -> Result<ClockModel, Report> {
      let record = runs.get(id)?;
      let file = record
        .output_files
        .iter()
        .find(|file| file.kind == OutputSelection::ClockModel)
        .expect("the run wrote a clock model");
      json_read_file(runs.store().out_dir(id).join(&file.path))
    }

    pub(super) fn finished_run(runs: &Arc<RunManager>, command: AppCommand, config: Value) -> JobId {
      let created = runs
        .create(CreateRunRequest {
          command,
          config,
          title: None,
          defer_start: true,
        })
        .unwrap();
      runs
        .start(&created.id, None, Box::new(|_config: &mut Value| Ok(())))
        .unwrap()
        .run();
      assert_eq!(RunStatus::Ok, runs.get(&created.id).unwrap().status);
      created.id
    }

    pub(super) fn zika(file: &str) -> PathBuf {
      dataset("zika/20", file)
    }

    pub(super) fn dataset(name: &str, file: &str) -> PathBuf {
      Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .join(name)
        .join(file)
    }

    pub(super) fn auspice_path(runs: &Arc<RunManager>, id: &JobId) -> PathBuf {
      let record = runs.get(id).unwrap();
      let file = record
        .output_files
        .iter()
        .find(|file| file.kind == OutputSelection::Auspice)
        .expect("the run wrote an Auspice file");
      runs.store().out_dir(id).join(&file.path)
    }
  }
}
