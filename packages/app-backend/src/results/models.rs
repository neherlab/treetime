use deser::Serialize;
use schemars::JsonSchema;
use std::path::PathBuf;
use treetime::gtr::get_gtr::{AlignmentModel, ModelEntry, ModelOutput, SequenceModel, TraitModel};
use treetime_schema::skip_serializing_optionals;

pub fn run_models(gtr: Option<&ModelOutput>) -> RunModels {
  gtr.map_or_else(RunModels::default, |gtr| RunModels {
    alignments: gtr.alignments.iter().map(AlignmentSummary::from_model).collect(),
    models: gtr
      .models
      .iter()
      .map(|entry| RunModel::from_entry(entry, &gtr.alignments))
      .collect(),
  })
}

/// The substitution models of a run, from its model document.
#[derive(Clone, Debug, Default, PartialEq, JsonSchema, Serialize)]
pub struct RunModels {
  /// Alignment files of the run, in the order of the `--alignment` list; empty for runs without alignments.
  pub alignments: Vec<AlignmentSummary>,
  /// Models in the order of the model document: one per partition, or one per trait.
  pub models: Vec<RunModel>,
}

/// One alignment file of a run.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct AlignmentSummary {
  /// Alignment number: the position of the file in the `--alignment` list, starting at 1.
  pub alignment: usize,
  /// Path of the alignment file.
  pub path: PathBuf,
  /// Number of sites of the alignment file.
  pub sites: usize,
  /// Rate of insertion and deletion events per unit of branch length in this file; absent when the run has no indel term.
  pub indel_rate: Option<f64>,
}

impl AlignmentSummary {
  fn from_model(model: &AlignmentModel) -> Self {
    Self {
      alignment: model.alignment,
      path: model.path.clone(),
      sites: model.sites,
      indel_rate: model.indel_rate,
    }
  }
}

/// A substitution model of a run, from its model document.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[schemars(tag = "kind", rename_all = "kebab-case")]
#[deser(tag = "kind", rename_all = "kebab-case")]
pub enum RunModel {
  /// Model of one partition: a set of sites of one alignment file.
  Sequence(SequenceModelSummary),
  /// Model of the states of one discrete trait.
  Trait(TraitModelSummary),
}

impl RunModel {
  fn from_entry(entry: &ModelEntry, alignments: &[AlignmentModel]) -> Self {
    match entry {
      ModelEntry::Sequence(model) => Self::Sequence(SequenceModelSummary::from_model(model, alignments)),
      ModelEntry::Trait(model) => Self::Trait(TraitModelSummary::from_model(model)),
    }
  }
}

/// The substitution model fitted to one partition.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct SequenceModelSummary {
  /// Partition number, starting at 1, in the order of the model document.
  pub partition: usize,
  /// Name of the partition entry; absent for partitions without a name.
  pub name: Option<String>,
  /// Number of the alignment file of the partition.
  pub alignment: usize,
  /// Path of the alignment file of the partition; absent when the model document does not list the file.
  pub path: Option<PathBuf>,
  /// Number of sites of the partition.
  pub sites: usize,
  /// Name of the model.
  pub model: String,
  /// Overall substitution rate.
  pub mu: f64,
  /// Characters of the states in the order of the rows and columns of `W` and the entries of `pi` in the model document.
  pub states: Vec<String>,
}

impl SequenceModelSummary {
  fn from_model(model: &SequenceModel, alignments: &[AlignmentModel]) -> Self {
    Self {
      partition: model.partition,
      name: model.name.clone(),
      alignment: model.alignment,
      path: alignments
        .iter()
        .find(|alignment| alignment.alignment == model.alignment)
        .map(|alignment| alignment.path.clone()),
      sites: model.sites,
      model: model.model.model_name().to_string(),
      mu: model.model.mu(),
      states: model.states.clone(),
    }
  }
}

/// The substitution model fitted to the states of one discrete trait.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
pub struct TraitModelSummary {
  /// Metadata column of the trait.
  pub attribute: String,
  /// Name of the model.
  pub model: String,
  /// Overall substitution rate.
  pub mu: f64,
  /// Trait states in the order of the rows and columns of `W` and the entries of `pi` in the model document.
  pub states: Vec<String>,
}

impl TraitModelSummary {
  fn from_model(model: &TraitModel) -> Self {
    Self {
      attribute: model.attribute.clone(),
      model: model.model.model_name().to_string(),
      mu: model.model.mu(),
      states: model.states.clone(),
    }
  }
}
