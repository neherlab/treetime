use crate::command::AppCommand;
use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
use crate::job::JobId;
use crate::results::clock::{ClockResults, clock_results};
use crate::results::methods::{Citation, citation, timetree_methods};
use crate::results::mugration::{MugrationResults, mugration_results};
use crate::results::mutations::{AncestralResults, ancestral_results};
use crate::results::outputs::{OutputProblem, RunOutputs};
use crate::results::timetree::{TimetreeOutputs, TimetreeResults, timetree_results};
use crate::results::tree::ResultTree;
use crate::runs::errors::conflict;
use crate::runs::manager::RunManager;
use crate::runs::record::{RunRecord, RunStatus};
use eyre::{Report, WrapErr};
use schemars::JsonSchema;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::Path;
use treetime::gtr::get_gtr::GtrOutput;
use treetime_utils::make_report;
use util_augur_node_data_json::AugurNodeDataJsonRefine;

/// Results of a finished run, read from its output files.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct RunResults {
  /// Tree of the run's Auspice file; absent when the run wrote none.
  pub tree: Option<ResultTree>,
  /// Results specific to the command of the run.
  pub results: CommandResults,
  /// Paragraph describing the analysis, for a methods section.
  pub methods: Option<String>,
  /// The publication to cite.
  pub citation: Citation,
  /// Output files that could not be read.
  pub problems: Vec<OutputProblem>,
}

/// Results specific to the command of a run.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "command", content = "data", rename_all = "kebab-case")]
pub enum CommandResults {
  Timetree(TimetreeResults),
  Clock(ClockResults),
  Ancestral(AncestralResults),
  Mugration(MugrationResults),
  Optimize(TreeSummary),
  Prune(TreeSummary),
}

/// Summary of a tree an `optimize` or `prune` run wrote.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct TreeSummary {
  /// Number of samples.
  pub samples: usize,
  /// Number of internal nodes.
  pub internal_nodes: usize,
  /// Number of nucleotide mutations on all branches.
  pub mutations: usize,
  /// Sum of the branch lengths in the node data file.
  pub total_branch_length: Option<f64>,
  /// Substitution model the run fitted.
  pub substitution_model: Option<SubstitutionModel>,
}

/// A fitted substitution model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SubstitutionModel {
  /// Name of the model.
  pub name: String,
  /// Overall substitution rate.
  pub mu: f64,
}

pub fn run_results(manager: &RunManager, id: &JobId) -> Result<RunResults, Report> {
  let record = finished_record(manager, id)?;
  results_of_record(&record, &manager.store().out_dir(id))
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

pub fn results_of_record(record: &RunRecord, out_dir: &Path) -> Result<RunResults, Report> {
  let RunOutputs {
    auspice,
    tree,
    clock_model,
    clock_rows,
    node_data,
    gtr,
    trace,
    coalescent,
    problems,
  } = RunOutputs::read(record.command, out_dir, &record.output_files);
  let (results, methods) = match record.command {
    AppCommand::Timetree => {
      let config: TreetimeTimetreeArgsRaw = run_config(record)?;
      let results = timetree_results(
        &TimetreeOutputs {
          tree: tree.as_ref(),
          clock_model: clock_model.as_ref(),
          clock_rows: clock_rows.as_deref(),
          node_data_clock: node_data.as_ref().and_then(|data| data.metadata.clock.as_ref()),
          trace: trace.as_deref().unwrap_or_default(),
          coalescent: coalescent.as_deref().unwrap_or_default(),
        },
        &config,
      );
      let methods = results
        .estimates
        .as_ref()
        .map(|estimates| timetree_methods(&record.treetime_version, &config, estimates));
      (CommandResults::Timetree(results), methods)
    },
    AppCommand::Clock => (
      CommandResults::Clock(clock_results(
        tree.as_ref(),
        clock_model.as_ref(),
        clock_rows.as_deref(),
      )),
      None,
    ),
    AppCommand::Ancestral => (CommandResults::Ancestral(ancestral_results(tree.as_ref())?), None),
    AppCommand::Mugration => {
      let config: TreetimeMugrationArgsRaw = run_config(record)?;
      let attribute = config
        .attribute
        .ok_or_else(|| make_report!("mugration run `{}` names no attribute", record.id.as_str()))?;
      (
        CommandResults::Mugration(mugration_results(auspice.as_ref(), tree.as_ref(), &attribute)?),
        None,
      )
    },
    AppCommand::Optimize => (
      CommandResults::Optimize(tree_summary(tree.as_ref(), node_data.as_ref(), gtr.as_ref())),
      None,
    ),
    AppCommand::Prune => (CommandResults::Prune(tree_summary(tree.as_ref(), None, None)), None),
  };
  Ok(RunResults {
    tree,
    results,
    methods,
    citation: citation(),
    problems,
  })
}

fn run_config<T: DeserializeOwned>(record: &RunRecord) -> Result<T, Report> {
  serde_json::from_value(Value::Object(record.config.clone()))
    .wrap_err_with(|| format!("When reading the configuration of run `{}`", record.id.as_str()))
}

fn tree_summary(
  tree: Option<&ResultTree>,
  node_data: Option<&AugurNodeDataJsonRefine>,
  gtr: Option<&GtrOutput>,
) -> TreeSummary {
  let samples = tree.map_or(0, |tree| tree.tips().count());
  TreeSummary {
    samples,
    internal_nodes: tree.map_or(0, |tree| tree.nodes.len() - samples),
    mutations: tree.map_or(0, |tree| tree.nodes.iter().map(|node| node.mutations.len()).sum()),
    total_branch_length: node_data
      .filter(|data| !data.nodes.is_empty())
      .map(|data| data.nodes.values().map(|node| node.branch_length).sum()),
    substitution_model: gtr.map(|gtr| SubstitutionModel {
      name: gtr.model_name().to_string(),
      mu: gtr.mu(),
    }),
  }
}
