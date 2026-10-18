use deser::Serialize;
use schemars::JsonSchema;
use std::path::PathBuf;
use treetime::gtr::get_gtr::{ModelEntry, ModelOutput};

pub fn run_models(gtr: Option<&ModelOutput>) -> Vec<RunModel> {
  gtr.map_or_else(Vec::new, |gtr| gtr.models.iter().map(RunModel::from_entry).collect())
}

/// A substitution model of a run, from its model document.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[schemars(tag = "kind", rename_all = "kebab-case")]
#[deser(tag = "kind", rename_all = "kebab-case")]
pub enum RunModel {
  /// Model of the sites of one alignment file.
  Sequence(SequenceModelSummary),
  /// Model of the states of one discrete trait.
  Trait(TraitModelSummary),
}

/// The substitution model fitted to one alignment file.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
pub struct SequenceModelSummary {
  /// Partition number: the position of the file in the `--alignment` list, starting at 1.
  pub partition: usize,
  /// Path of the alignment file.
  pub alignment: PathBuf,
  /// Name of the model.
  pub name: String,
  /// Overall substitution rate.
  pub mu: f64,
}

/// The substitution model fitted to the states of one discrete trait.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
pub struct TraitModelSummary {
  /// Metadata column of the trait.
  pub attribute: String,
  /// Name of the model.
  pub name: String,
  /// Overall substitution rate.
  pub mu: f64,
  /// Trait states in the order of the rows and columns of `W` and the entries of `pi` in the model document.
  pub states: Vec<String>,
}

impl RunModel {
  fn from_entry(entry: &ModelEntry) -> Self {
    match entry {
      ModelEntry::Sequence(model) => Self::Sequence(SequenceModelSummary {
        partition: model.partition,
        alignment: model.alignment.clone(),
        name: model.model.model_name().to_string(),
        mu: model.model.mu(),
      }),
      ModelEntry::Trait(model) => Self::Trait(TraitModelSummary {
        attribute: model.attribute.clone(),
        name: model.model.model_name().to_string(),
        mu: model.model.mu(),
        states: model.states.clone(),
      }),
    }
  }
}
