use crate::command::{OutputAlignment, OutputFile};
use crate::job::JobId;
use crate::results::citation::{Citation, citation};
use crate::results::clock::{ClockResults, clock_results};
use crate::results::homoplasy::{HomoplasyResults, homoplasy_results};
use crate::results::models::{RunModels, run_models};
use crate::results::mugration::{MugrationResults, mugration_results};
use crate::results::mutations::{AncestralResults, ancestral_results};
use crate::results::outputs::{OutputProblem, RunOutputs, shown_auspice};
use crate::results::timetree::{TimetreeOutputs, TimetreeResults, timetree_results};
use crate::results::tree::ResultTree;
use crate::results::tree_change::{TreeChange, TreeStatistics, input_tree_statistics, tree_change};
use crate::runs::errors::conflict;
use crate::runs::input_store::InputStore;
use crate::runs::manager::RunManager;
use crate::runs::record::{RunRecord, RunStatus};
use crate::runs::run_paths::{OUT_DIR, resolve_path};
use app_commands::command_config::CommandConfig;
use app_commands::commands::shared::tree_input::{TreeArgs, TreeDialectArgs};
use deser::Serialize;
use eyre::Report;
use schemars::JsonSchema;
use std::path::Path;
use treetime::alphabet::alphabet::Alphabet;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::make_report;

pub fn run_results(manager: &RunManager, id: &JobId) -> Result<RunResults, Report> {
  let record = finished_record(manager, id)?;
  results_of_record(&record, &manager.store().run_dir(id), manager.input_store())
}

pub fn finished_record(manager: &RunManager, id: &JobId) -> Result<RunRecord, Report> {
  let record = manager.get(id)?;
  if record.status != RunStatus::Ok {
    return Err(conflict(format!(
      "run `{}` has no results: it is {}",
      id.as_str(),
      record.status
    )));
  }
  Ok(record)
}

pub fn results_of_record(record: &RunRecord, run_dir: &Path, store: &InputStore) -> Result<RunResults, Report> {
  let out_dir = run_dir.join(OUT_DIR);
  let RunOutputs {
    auspice,
    tree,
    clock_model,
    clock_rows,
    node_data,
    gtr,
    trace,
    coalescent,
    homoplasy,
    mut problems,
  } = RunOutputs::read(record.config.command(), &out_dir, &record.output_files);
  let input = InputTree { record, run_dir, store };
  let results = match &record.config {
    CommandConfig::Timetree { config } => CommandResults::Timetree(Box::new(timetree_results(
      &TimetreeOutputs {
        tree: tree.as_ref(),
        clock_model: clock_model.as_ref(),
        clock_rows: clock_rows.as_deref(),
        node_data_clock: node_data.as_ref().and_then(|data| data.metadata.clock.as_ref()),
        trace: trace.as_deref().unwrap_or_default(),
        coalescent: coalescent.as_deref().unwrap_or_default(),
      },
      config,
    ))),
    CommandConfig::Clock { .. } => CommandResults::Clock(clock_results(
      tree.as_ref(),
      clock_model.as_ref(),
      clock_rows.as_deref(),
    )),
    CommandConfig::Ancestral { config } => CommandResults::Ancestral(ancestral_results(
      tree.as_ref(),
      &Alphabet::new(config.alphabet_args.alphabet_name().unwrap_or_default())?,
    )?),
    CommandConfig::Homoplasy { .. } => CommandResults::Homoplasy(Box::new(homoplasy_results(homoplasy.as_ref()))),
    CommandConfig::Mugration { config } => {
      let attribute = config
        .attribute
        .clone()
        .ok_or_else(|| make_report!("mugration run `{}` names no attribute", record.id.as_str()))?;
      CommandResults::Mugration(mugration_results(auspice.as_ref(), tree.as_ref(), &attribute)?)
    },
    CommandConfig::Optimize { config } => CommandResults::Optimize(tree_change(
      input.statistics(&config.tree, config.newick_dialect, &mut problems),
      tree.as_ref(),
    )),
    CommandConfig::Prune { config } => CommandResults::Prune(tree_change(
      input.statistics(&config.tree, config.newick_dialect, &mut problems),
      tree.as_ref(),
    )),
  };
  Ok(RunResults {
    tree,
    tree_alignment: tree_alignment(&record.output_files),
    models: run_models(gtr.as_ref()),
    results,
    citation: citation(),
    problems,
  })
}

/// Results of a finished run, read from its output files.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct RunResults {
  /// Tree of the run's Auspice file; absent when the run wrote none.
  pub tree: Option<ResultTree>,
  /// Alignment file whose Auspice file gives the tree, when the run has several alignment files.
  pub tree_alignment: Option<OutputAlignment>,
  /// Substitution models and alignment files of the run's model document; empty when the run wrote none.
  pub models: RunModels,
  /// Results specific to the command of the run.
  pub results: CommandResults,
  /// The publication to cite.
  pub citation: Citation,
  /// Output files, and the input tree of an `optimize` or `prune` run, that could not be read.
  pub problems: Vec<OutputProblem>,
}

/// Results specific to the command of a run.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[schemars(tag = "command", content = "data", rename_all = "kebab-case")]
#[deser(tag = "command", content = "data", rename_all = "kebab-case")]
pub enum CommandResults {
  Timetree(Box<TimetreeResults>),
  Clock(ClockResults),
  Ancestral(AncestralResults),
  Homoplasy(Box<HomoplasyResults>),
  Mugration(MugrationResults),
  Optimize(TreeChange),
  Prune(TreeChange),
}

fn tree_alignment(output_files: &[OutputFile]) -> Option<OutputAlignment> {
  let several = output_files
    .iter()
    .any(|file| file.alignment.as_ref().is_some_and(|alignment| alignment.number > 1));
  shown_auspice(output_files)
    .filter(|_| several)
    .and_then(|file| file.alignment.clone())
}

struct InputTree<'a> {
  record: &'a RunRecord,
  run_dir: &'a Path,
  store: &'a InputStore,
}

impl InputTree<'_> {
  fn statistics(
    &self,
    tree: &TreeArgs,
    dialect: TreeDialectArgs,
    problems: &mut Vec<OutputProblem>,
  ) -> Option<TreeStatistics> {
    let path = resolve_path(tree.tree.as_deref()?, &self.record.inputs, self.run_dir, self.store);
    match input_tree_statistics(&path, dialect.newick_dialect.dialect()) {
      Ok(statistics) => Some(statistics),
      Err(report) => {
        problems.push(OutputProblem {
          path: path.to_string_lossy().into_owned(),
          message: format!("{report:#}"),
        });
        None
      },
    }
  }
}
