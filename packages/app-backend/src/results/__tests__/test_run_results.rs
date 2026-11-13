#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::sparse;
  use crate::command::OutputAlignment;
  use crate::config::catalog::SettingPurpose;
  use crate::json_float::JsonFloat;
  use crate::results::auspice::run_auspice;
  use crate::results::clades::{CladeRequest, clade_in_runs};
  use crate::results::clock::{ClockLine, RootToTip};
  use crate::results::compare::{EstimateDelta, LikelihoodMismatch, compare_runs, likelihood_mismatch_text};
  use crate::results::models::{RunModel, RunModels};
  use crate::results::outputs::OutputProblem;
  use crate::results::run_results::{CommandResults, RunResults, run_results};
  use crate::results::tree::ResultTree;
  use crate::runs::headline::RunHeadline;
  use crate::runs::manager::RunManager;
  use crate::runs::record::CreateRunRequest;
  use crate::runs::setting_differences::{DifferenceValues, SettingDifference};
  use app_commands::command::AppCommand;
  use app_commands::json_value::JsonValue;
  use app_output::output_plan::OutputSelection;
  use eyre::Report;
  use helpers::{auspice_path, clock_model, finished_run, output_path, timetree_config, zika};
  use itertools::{Itertools, izip};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::fs;
  use tempfile::tempdir;
  use treetime::clock::rtt::ClockDateSource;
  use treetime::o;
  use treetime_utils::datetime::year_fraction::year_fraction_days_between;
  use treetime_utils::io::json::{json_read_file, to_json_value};
  use treetime_utils::{assert_error, pretty_assert_relative_eq, pretty_assert_ulps_eq, vec_of_owned};

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
  fn test_run_results_headline_of_a_clock_run_holds_the_root_date_of_the_regression() -> Result<(), Report> {
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
    let estimates = clock.estimates;
    assert_eq!(
      (
        true,
        RunHeadline {
          root_date: estimates.root_date.clone(),
          clock_rate: estimates.clock_rate.map(JsonFloat),
          r_squared: estimates.r_squared.map(JsonFloat),
          recurrent_substitutions: None,
        }
      ),
      (estimates.root_date.is_some(), runs.get(&id)?.headline)
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
        "pi": [0.5, 0.5],
        "rates": [[0.0, mu], [mu, 0.0]],
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
  fn test_run_results_of_a_prune_run_describe_the_tree_before_and_after() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let inputs = tempdir()?;
    let tree = inputs.path().join("tree.nwk");
    fs::write(&tree, "((A:0.1,B:0.1):0.1,(C:0.2,D:0.3):0.05,E:0.3);\n")?;
    let id = finished_run(&runs, AppCommand::Prune, json!({ "tree": tree, "remove": ["E"] }));

    let CommandResults::Prune(change) = run_results(&runs, &id)?.results else {
      panic!("a prune run has prune results");
    };

    let [samples, internal, polytomies, length] = change.rows.as_slice() else {
      panic!("a tree change has four statistics");
    };
    assert_eq!(
      vec![
        (o!("Samples"), true, Some(5.0), Some(4.0)),
        (o!("Internal nodes"), true, Some(3.0), Some(3.0)),
        (o!("Polytomies"), true, Some(1.0), Some(0.0)),
      ],
      [samples, internal, polytomies]
        .map(|row| (row.label.clone(), row.count, row.before, row.after))
        .to_vec()
    );
    assert_eq!((o!("Total branch length"), false), (length.label.clone(), length.count));
    pretty_assert_ulps_eq!(
      1.15,
      length.before.expect("the input tree has branch lengths"),
      max_ulps = 4
    );
    pretty_assert_ulps_eq!(
      0.85,
      length.after.expect("the result tree has branch lengths"),
      max_ulps = 4
    );
    Ok(())
  }

  #[test]
  fn test_run_results_of_a_prune_run_name_an_input_tree_that_is_gone() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let inputs = tempdir()?;
    let tree = inputs.path().join("tree.nwk");
    fs::write(&tree, "((A:0.1,B:0.1):0.1,C:0.3);\n")?;
    let id = finished_run(&runs, AppCommand::Prune, json!({ "tree": tree, "remove": ["C"] }));
    fs::remove_file(&tree)?;

    let results = run_results(&runs, &id)?;

    let CommandResults::Prune(change) = &results.results else {
      panic!("a prune run has prune results");
    };
    assert_eq!(
      (
        None,
        vec![OutputProblem {
          path: tree.display().to_string(),
          message: format!(
            "When opening file '{}': No such file or directory (os error 2)",
            tree.display()
          ),
        }]
      ),
      (change.rows.iter().find_map(|row| row.before), results.problems)
    );
    Ok(())
  }

  #[test]
  fn test_run_files_list_the_main_formats_of_a_time_tree_first_in_menu_order() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let id = finished_run(&runs, AppCommand::Timetree, timetree_config(None));

    let files = runs.files(&id)?;

    let main = files
      .iter()
      .take_while(|file| file.main)
      .map(|file| (file.path.as_str(), file.label.as_str()))
      .collect::<Vec<_>>();
    assert_eq!(
      (
        vec![
          ("timetree.auspice.json", "Auspice JSON"),
          ("timetree.nwk", "Dated tree (Newick, branch lengths in years)"),
          ("timetree.nexus", "Dated tree with node annotations (Nexus)"),
        ],
        true
      ),
      (main.clone(), files.iter().skip(main.len()).all(|file| !file.main))
    );
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
      source_run_id: None,
      example: None,
    })?;

    assert_error!(
      run_results(&runs, &created.id),
      format!("run `{}` has no results: it is created", created.id.as_str())
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
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
      source_run_id: None,
      example: None,
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
    assert_eq!(
      Some(EstimateDelta {
        root_date_days: Some(0.0),
        clock_rate_percent: Some(0.0),
        summary: o!("Root date 0 d, clock rate 0.0 %"),
      }),
      comparison.delta
    );
    let estimates = comparison.estimates.expect("both runs are time trees");
    let ancestors = comparison.ancestors.expect("both runs are time trees");
    assert_eq!(
      (true, true, true),
      (
        comparison.same_analysis,
        settings.differences.is_empty(),
        settings.same_config
      )
    );
    assert_eq!(
      (Some(0.0), Some(0.0), 0, Some(0.0), vec![], ancestors.ancestors, true),
      (
        estimates.root_shift_days,
        estimates.clock_rate_change_percent,
        estimates.excluded_samples_change,
        estimates.log_likelihood_change,
        estimates.log_likelihood_mismatches,
        ancestors.shifts.len(),
        ancestors.shifts.iter().all(|shift| shift.shift_days == 0.0)
      )
    );
    Ok(())
  }

  #[test]
  fn test_compare_runs_with_different_drawn_seeds_have_different_configs() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let mut config = timetree_config(None);
    config
      .as_object_mut()
      .expect("a configuration is a mapping")
      .remove("seed");
    let first = finished_run(&runs, AppCommand::Timetree, config.clone());
    let second = finished_run(&runs, AppCommand::Timetree, config);

    let settings = compare_runs(&runs, &first, &second)?
      .settings
      .expect("both runs execute the same command");

    assert_eq!(
      (vec![o!("seed")], false),
      (
        settings
          .differences
          .iter()
          .map(|difference| difference.key.clone())
          .collect_vec(),
        settings.same_config
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
        vec![SettingDifference {
          key: o!("clock_rate"),
          label: o!("Clock rate"),
          group: o!("Molecular clock"),
          purpose: SettingPurpose::Analysis,
          values: DifferenceValues::Setting {
            first: None,
            second: Some(JsonValue(json!(1e-3))),
          },
        }],
        false,
        true,
      ),
      (settings.differences, settings.same_config, settings.compared > 1)
    );
    Ok(())
  }

  #[test]
  fn test_compare_runs_with_different_coalescent_priors_leaves_out_the_log_likelihood_change() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let first = finished_run(&runs, AppCommand::Timetree, timetree_config(None));
    let mut optimized = timetree_config(None);
    optimized["coalescent"] = json!("opt");
    let second = finished_run(&runs, AppCommand::Timetree, optimized);

    let estimates = compare_runs(&runs, &first, &second)?
      .estimates
      .expect("both runs are time trees");

    assert_eq!(
      (true, None),
      (
        estimates
          .log_likelihood_mismatches
          .contains(&LikelihoodMismatch::CoalescentPrior),
        estimates.log_likelihood_change
      )
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::comparable(vec![],                                                                None)]
  #[case::one(       vec![LikelihoodMismatch::CoalescentPrior],                             Some("different coalescent prior"))]
  #[case::two(       vec![LikelihoodMismatch::Data, LikelihoodMismatch::Topology],          Some("different data, different tree topology"))]
  #[trace]
  fn test_compare_likelihood_mismatch_text_names_each_mismatch(
    #[case] mismatches: Vec<LikelihoodMismatch>,
    #[case] expected: Option<&str>,
  ) {
    assert_eq!(expected.map(str::to_owned), likelihood_mismatch_text(&mismatches));
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
      (false, None, None, None, None),
      (
        comparison.same_analysis,
        comparison.settings,
        comparison.delta,
        comparison.estimates,
        comparison.ancestors
      )
    );
    Ok(())
  }

  #[test]
  fn test_compare_clock_runs_gives_the_change_of_the_regression_root_date_and_clock_rate() -> Result<(), Report> {
    let root = tempdir()?;
    let runs = RunManager::open(root.path())?;
    let config =
      |keep_root: bool| json!({ "tree": zika("tree.nwk"), "metadata": zika("metadata.tsv"), "keep_root": keep_root });
    let first = finished_run(&runs, AppCommand::Clock, config(false));
    let second = finished_run(&runs, AppCommand::Clock, config(true));
    let (first_model, second_model) = (clock_model(&runs, &first)?, clock_model(&runs, &second)?);
    let days =
      year_fraction_days_between(first_model.t_mrca(), second_model.t_mrca()).expect("both dates are in range");
    let percent = (second_model.clock_rate() / first_model.clock_rate() - 1.0) * 100.0;

    let delta = compare_runs(&runs, &first, &second)?
      .delta
      .expect("both runs are clock runs");

    assert_eq!(
      (format!("Root date {days:+.0} d, clock rate {percent:+.1} %"), true),
      (delta.summary.clone(), days.round() != 0.0)
    );
    pretty_assert_ulps_eq!(
      days,
      delta.root_date_days.expect("both runs have a root date"),
      max_ulps = 4
    );
    pretty_assert_relative_eq!(
      percent,
      delta.clock_rate_percent.expect("both runs have a clock rate"),
      max_relative = 1e-9
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
          source_run_id: None,
          example: None,
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
