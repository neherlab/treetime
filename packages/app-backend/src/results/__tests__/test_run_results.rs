#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::sparse;
  use crate::command::OutputAlignment;
  use crate::json_float::JsonFloat;
  use crate::results::auspice::run_auspice;
  use crate::results::clades::{CladeRequest, clade_in_runs};
  use crate::results::clock::{ClockLine, RootToTip};
  use crate::results::compare::compare_runs;
  use crate::results::models::{RunModel, RunModels};
  use crate::results::run_results::{CommandResults, RunResults, run_results};
  use crate::results::tree::ResultTree;
  use crate::runs::headline::RunHeadline;
  use crate::runs::manager::RunManager;
  use crate::runs::record::CreateRunRequest;
  use crate::runs::setting_differences::SettingDifference;
  use app_commands::command::AppCommand;
  use app_commands::json_value::JsonValue;
  use app_output::output_plan::OutputSelection;
  use eyre::Report;
  use helpers::{
    auspice_path, clock_model, dataset, finished_run, homoplasy_config, output_path, timetree_config, zika,
  };
  use itertools::izip;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::collections::BTreeSet;
  use std::fs;
  use tempfile::tempdir;
  use treetime::clock::rtt::ClockDateSource;
  use treetime::o;
  use treetime_utils::io::json::{json_read_file, to_json_value};
  use treetime_utils::{assert_error, vec_of_owned};

  #[test]
  fn test_run_results_show_the_tree_of_the_first_alignment_file() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let inputs = tempdir()?;
    let tree = inputs.path().join("tree.nwk");
    let first = inputs.path().join("ha.fasta");
    let second = inputs.path().join("na.fasta");
    fs::write(&tree, "((A:0.1,B:0.1):0.1,C:0.1);\n")?;
    fs::write(&first, ">A\nACGT\n>B\nACGT\n>C\nGCGT\n")?;
    fs::write(&second, ">A\nTTAA\n>B\nTAAA\n>C\nTTAA\n")?;
    let id = finished_run(
      &runs,
      AppCommand::Ancestral,
      json!({ "tree": tree, "alignment": [&second, &first] }),
    );

    let record = runs.get(&id)?;
    let out_dir = runs.store().out_dir(&id);
    let auspice_of = |number: usize| -> Result<ResultTree, Report> {
      let file = record
        .output_files
        .iter()
        .find(|file| {
          file.kind == OutputSelection::Auspice
            && file
              .alignment
              .as_ref()
              .is_some_and(|alignment| alignment.number == number)
        })
        .expect("the run wrote one Auspice file per alignment file");
      ResultTree::from_auspice(&json_read_file(out_dir.join(&file.path))?)
    };
    let results = run_results(&runs, &id)?;
    let shown = ResultTree::from_auspice(&run_auspice(&runs, &id)?.0)?;
    let listed_auspice = runs
      .files(&id)?
      .into_iter()
      .filter(|file| file.kind == Some(OutputSelection::Auspice))
      .map(|file| (file.path, file.alignment))
      .collect::<Vec<_>>();

    assert_eq!(
      vec![
        (
          o!("ancestral.1.auspice.json"),
          Some(OutputAlignment {
            number: 1,
            path: second.clone()
          })
        ),
        (
          o!("ancestral.2.auspice.json"),
          Some(OutputAlignment { number: 2, path: first })
        ),
      ],
      listed_auspice
    );
    assert_eq!(
      (
        true,
        Some(auspice_of(1)?),
        auspice_of(1)?,
        Some(OutputAlignment {
          number: 1,
          path: second
        }),
      ),
      (
        auspice_of(1)? != auspice_of(2)?,
        results.tree,
        shown,
        results.tree_alignment,
      )
    );
    Ok(())
  }

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
        recurrent_substitutions: None,
      },
      runs.get(&id)?.headline
    );
    Ok(())
  }

  #[test]
  fn test_run_results_of_a_homoplasy_run_read_its_statistics_file() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, AppCommand::Homoplasy, homoplasy_config());

    let RunResults {
      tree,
      results: CommandResults::Homoplasy(homoplasy),
      problems,
      ..
    } = run_results(&runs, &id)?
    else {
      panic!("a homoplasy run has homoplasy results");
    };

    let record = runs.get(&id)?;
    let homoplasy = homoplasy.statistics.expect("the run wrote its statistics file");
    let kinds = record
      .output_files
      .iter()
      .map(|file| file.kind)
      .collect::<BTreeSet<_>>();
    assert_eq!(
      (
        (true, true),
        Vec::<String>::new(),
        20,
        homoplasy.sites_hit_more_than_once,
        Some(homoplasy.recurrent_substitutions),
      ),
      (
        (
          kinds.contains(&OutputSelection::Auspice),
          kinds.contains(&OutputSelection::HomoplasyStats)
        ),
        problems.into_iter().map(|problem| problem.message).collect(),
        tree.expect("the run wrote an Auspice tree").tips().count(),
        homoplasy.sites.len(),
        record.headline.recurrent_substitutions,
      )
    );
    Ok(())
  }

  #[test]
  fn test_run_results_of_a_homoplasy_run_without_statistics_file_are_absent() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, AppCommand::Homoplasy, homoplasy_config());
    let record = runs.get(&id)?;
    let stats = record
      .output_files
      .iter()
      .find(|file| file.kind == OutputSelection::HomoplasyStats)
      .expect("the run wrote its statistics file");
    fs::write(runs.store().out_dir(&id).join(&stats.path), "{}")?;

    let RunResults {
      results: CommandResults::Homoplasy(homoplasy),
      problems,
      ..
    } = run_results(&runs, &id)?
    else {
      panic!("a homoplasy run has homoplasy results");
    };

    assert_eq!(
      (None, vec![stats.path.to_string_lossy().into_owned()]),
      (
        homoplasy.statistics,
        problems.into_iter().map(|problem| problem.path).collect::<Vec<_>>()
      )
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
  fn test_run_results_of_an_optimize_run_list_one_substitution_model_per_partition_and_the_indel_rate_per_alignment()
  -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(
      &runs,
      AppCommand::Optimize,
      json!({ "tree": zika("tree.nwk"), "alignment": [zika("aln.fasta.xz")] }),
    );
    let model = |mu: f64| {
      json!({
        "model_type": "inferred",
        "model_name": "infer",
        "mu": mu,
        "pi": [0.5, 0.5],
        "W": [[0.0, 2.0], [2.0, 0.0]],
        "n_states": 2,
      })
    };
    let entry = |partition: usize, name: &str, alignment: usize, mu: f64| {
      let mut entry = json!({
        "kind": "sequence",
        "partition": partition,
        "name": name,
        "alignment": alignment,
        "sites": 50,
        "states": ["A", "T"],
      });
      entry
        .as_object_mut()
        .expect("entry is an object")
        .extend(model(mu).as_object().expect("model is an object").clone());
      entry
    };
    fs::write(
      output_path(&runs, &id, OutputSelection::Model),
      json!({
        "alignments": [
          { "alignment": 1, "path": "ha/aligned.fasta", "sites": 100, "indel_rate": 0.5 },
          { "alignment": 2, "path": "na/aligned.fasta", "sites": 50 },
        ],
        "models": [entry(1, "cp12", 1, 0.75), entry(2, "cp3", 1, 1.25), entry(3, "na", 2, 1.0)],
      })
      .to_string(),
    )?;

    let results = run_results(&runs, &id)?;

    let summary = |partition: usize, name: &str, alignment: usize, path: &str, mu: f64| {
      json!({
        "kind": "sequence",
        "partition": partition,
        "name": name,
        "alignment": alignment,
        "path": path,
        "sites": 50,
        "model": "Infer",
        "mu": mu,
        "states": ["A", "T"],
      })
    };
    assert_eq!(
      json!({
        "alignments": [
          { "alignment": 1, "path": "ha/aligned.fasta", "sites": 100, "indel_rate": 0.5 },
          { "alignment": 2, "path": "na/aligned.fasta", "sites": 50 },
        ],
        "models": [
          summary(1, "cp12", 1, "ha/aligned.fasta", 0.75),
          summary(2, "cp3", 1, "ha/aligned.fasta", 1.25),
          summary(3, "na", 2, "na/aligned.fasta", 1.0),
        ],
      }),
      to_json_value(&results.models)?
    );
    Ok(())
  }

  #[rstest]
  #[case::ancestral_marginal(AppCommand::Ancestral, json!({ "tree": zika("tree.nwk"), "alignment": [zika("aln.fasta.xz")], "method_anc": "marginal" }))]
  fn test_run_results_list_the_sequence_model_the_run_fitted(
    #[case] command: AppCommand,
    #[case] config: Value,
  ) -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, command, config);

    let models = run_results(&runs, &id)?.models.models;

    let [RunModel::Sequence(model)] = models.as_slice() else {
      panic!("the run lists one sequence model, got {models:?}");
    };
    assert_eq!(
      (
        1,
        1,
        Some(zika("aln.fasta.xz")),
        vec_of_owned!["A", "C", "G", "T"],
        true
      ),
      (
        model.partition,
        model.alignment,
        model.path.clone(),
        model.states.clone(),
        model.mu > 0.0
      )
    );
    Ok(())
  }

  #[test]
  fn test_run_results_of_a_mugration_run_list_the_trait_model_with_its_states() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(
      &runs,
      AppCommand::Mugration,
      json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv"), "attribute": "country" }),
    );

    let results = run_results(&runs, &id)?;

    let CommandResults::Mugration(mugration) = &results.results else {
      panic!("a mugration run has mugration results");
    };
    let [RunModel::Trait(model)] = results.models.models.as_slice() else {
      panic!("a mugration run lists one trait model, got {:?}", results.models.models);
    };
    assert_eq!(
      ("country", true, true),
      (
        model.attribute.as_str(),
        model.states.len() >= mugration.states,
        model.mu > 0.0
      )
    );
    Ok(())
  }

  #[test]
  fn test_run_results_of_a_fitch_ancestral_run_list_no_model() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(
      &runs,
      AppCommand::Ancestral,
      json!({ "tree": zika("tree.nwk"), "alignment": [zika("aln.fasta.xz")], "method_anc": "parsimony" }),
    );

    let results = run_results(&runs, &id)?;

    assert_eq!((RunModels::default(), vec![]), (results.models, results.problems));
    Ok(())
  }

  #[test]
  fn test_run_results_refuse_a_run_that_did_not_finish() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let created = runs.create(&CreateRunRequest {
      command: AppCommand::Clock,
      config: sparse(json!({ "tree": zika("tree.nwk") })),
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

    let document = to_json_value(&run_auspice(&runs, &id)?)?;

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
    let expected = izip!(&region.states, ["#332288", "#88ccee", "#44aa99", "#117733"])
      .map(|(state, color)| [state.as_str(), color])
      .collect::<Vec<_>>();
    assert_eq!(Some(json!(expected)), scale);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::mugration(AppCommand::Mugration, json!({ "tree": dataset("dengue/20", "tree.nwk"), "metadata": dataset("dengue/20", "metadata.tsv"), "metadata_id_columns": ["genbank_accession"], "attribute": "region" }))]
  #[case::timetree( AppCommand::Timetree,  timetree_config(None))]
  #[case::ancestral(AppCommand::Ancestral, json!({ "tree": zika("tree.nwk"), "alignment": [zika("aln.fasta.xz")] }))]
  #[case::clock(    AppCommand::Clock,     json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv") }))]
  #[case::homoplasy(AppCommand::Homoplasy, homoplasy_config())]
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

    let document = to_json_value(&run_auspice(&runs, &id)?)?;

    let mut expected: Value = serde_json::from_slice(&written)?;
    let displayed = document["meta"]["colorings"].as_array().expect("colorings");
    for (coloring, shown) in izip!(expected["meta"]["colorings"].as_array_mut().expect("colorings"), displayed) {
      if let Some(scale) = shown.get("scale") {
        coloring["scale"] = scale.clone();
      }
    }
    assert_eq!((written, expected), (fs::read(&path)?, document));
    Ok(())
  }

  #[test]
  fn test_run_auspice_refuses_a_run_that_did_not_finish() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let created = runs.create(&CreateRunRequest {
      command: AppCommand::Clock,
      config: sparse(json!({ "tree": zika("tree.nwk") })),
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
          first: None,
          second: Some(JsonValue(json!(1e-3))),
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
    use crate::__tests__::test_support::tests::sparse;
    use crate::job::JobId;
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, RunStatus, StartRunRequest};
    use app_commands::command::AppCommand;
    use app_output::output_plan::OutputSelection;
    use eyre::Report;
    use serde_json::{Value, json};
    use std::path::{Path, PathBuf};
    use std::sync::Arc;
    use treetime::clock::clock_model::ClockModel;
    use treetime_utils::io::json::json_read_file;

    pub(super) fn timetree_config(clock_rate: Option<f64>) -> Value {
      let mut config = json!({
        "tree": zika("tree.nwk"),
        "metadata": zika("metadata.tsv"),
        "alignment": [zika("aln.fasta.xz")],
        "max_iter": 2,
        "seed": 7,
      });
      if let Some(clock_rate) = clock_rate {
        config["clock_rate"] = json!(clock_rate);
      }
      config
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
        .create(&CreateRunRequest {
          command,
          config: sparse(config),
          defer_start: true,
        })
        .unwrap();
      runs
        .start(
          &created.id,
          &StartRunRequest::default(),
          Box::new(|_config: &mut Value| Ok(())),
        )
        .unwrap()
        .run();
      assert_eq!(RunStatus::Ok, runs.get(&created.id).unwrap().status);
      created.id
    }

    pub(super) fn homoplasy_config() -> Value {
      json!({ "tree": zika("tree.nwk"), "alignment": [zika("aln.fasta.xz")] })
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
      output_path(runs, id, OutputSelection::Auspice)
    }

    pub(super) fn output_path(runs: &Arc<RunManager>, id: &JobId, kind: OutputSelection) -> PathBuf {
      let record = runs.get(id).unwrap();
      let file = record
        .output_files
        .iter()
        .find(|file| file.kind == kind)
        .expect("the run wrote the output file");
      runs.store().out_dir(id).join(&file.path)
    }
  }
}
