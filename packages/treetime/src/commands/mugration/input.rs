use crate::partition::payload::ancestral::GraphAncestral;
use std::collections::BTreeMap;

/// Parsed mugration input for the execution layer.
///
/// Represents already-parsed tree and trait data plus mugration configuration,
/// without embedding output destinations or file paths.
#[derive(Debug)]
pub struct MugrationInput {
  /// In-memory tree.
  pub graph: GraphAncestral,

  /// Trait values keyed by leaf name.
  pub traits: BTreeMap<String, String>,

  /// Attribute name being reconstructed (e.g., "country").
  pub attribute: String,

  /// Optional weights map keyed by state name.
  pub weights: Option<BTreeMap<String, f64>>,

  /// String indicating missing data (e.g., "?").
  pub missing_data: String,

  /// Optional pseudo-counts for equilibrium frequency smoothing.
  pub pc: Option<f64>,

  /// Threshold for portion of attribute values allowed to be missing from weights.
  pub missing_weights_threshold: f64,

  /// Number of GTR refinement iterations (re-estimation of rate matrix from data).
  pub iterations: usize,

  /// Factor to inflate overall switching rate to counteract sampling bias.
  pub sampling_bias_correction: Option<f64>,
}

/// Parameters for mugration model construction.
///
/// Extracted from MugrationInput for internal use during model setup.
#[derive(Clone, Debug)]
pub struct MugrationParams {
  /// Attribute name being reconstructed.
  pub attribute: String,

  /// String indicating missing data.
  pub missing_data: String,

  /// Optional pseudo-counts.
  pub pc: Option<f64>,

  /// Threshold for missing weights validation.
  pub missing_weights_threshold: f64,

  /// Number of GTR refinement iterations.
  pub iterations: usize,

  /// Sampling bias correction factor.
  pub sampling_bias_correction: Option<f64>,
}

impl MugrationInput {
  pub fn params(&self) -> MugrationParams {
    MugrationParams {
      attribute: self.attribute.clone(),
      missing_data: self.missing_data.clone(),
      pc: self.pc,
      missing_weights_threshold: self.missing_weights_threshold,
      iterations: self.iterations,
      sampling_bias_correction: self.sampling_bias_correction,
    }
  }
}
