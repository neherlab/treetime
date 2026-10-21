use crate::alphabet::alphabet::Alphabet;
use crate::constants::POSTERIOR_TIE_TOLERANCE;
use derive_more::{Display, Error};
use deser::adapters::TryFromInto;
use deser::{Deserialize, Serialize};
use ndarray::{ArrayBase, Data, Ix1};
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use std::borrow::Cow;
use std::str::FromStr;
use treetime_primitives::{AsciiChar, StateSet};
use treetime_utils::array::ndarray::{argmax_all_within, argmax_first_within};

pub fn most_likely_state<S: Data<Elem = f64>>(posterior: &ArrayBase<S, Ix1>) -> Option<usize> {
  argmax_first_within(posterior, POSTERIOR_TIE_TOLERANCE)
}

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
#[deser(as = TryFromInto<f64>)]
pub struct AmbiguityThreshold(f64);

impl AmbiguityThreshold {
  pub fn new(value: f64) -> Result<Self, InvalidAmbiguityThreshold> {
    if (0.0..=1.0).contains(&value) {
      Ok(Self(value))
    } else {
      Err(InvalidAmbiguityThreshold { value })
    }
  }

  pub fn get(self) -> f64 {
    self.0
  }

  pub fn plausible_states<S: Data<Elem = f64>>(
    self,
    posterior: &ArrayBase<S, Ix1>,
  ) -> impl Iterator<Item = usize> + '_ {
    argmax_all_within(posterior, self.0 + POSTERIOR_TIE_TOLERANCE)
  }

  pub(crate) fn ambiguity_char<S: Data<Elem = f64>>(
    self,
    posterior: &ArrayBase<S, Ix1>,
    alphabet: &Alphabet,
  ) -> Option<AsciiChar> {
    let states: StateSet = self
      .plausible_states(posterior)
      .map(|index| alphabet.char(index))
      .collect();
    (states.len() > 1).then(|| alphabet.ambiguity_char(states))
  }
}

impl TryFrom<f64> for AmbiguityThreshold {
  type Error = InvalidAmbiguityThreshold;

  fn try_from(value: f64) -> Result<Self, InvalidAmbiguityThreshold> {
    Self::new(value)
  }
}

impl From<AmbiguityThreshold> for f64 {
  fn from(threshold: AmbiguityThreshold) -> Self {
    threshold.0
  }
}

impl FromStr for AmbiguityThreshold {
  type Err = InvalidAmbiguityThresholdText;

  fn from_str(text: &str) -> Result<Self, InvalidAmbiguityThresholdText> {
    let value = text
      .trim()
      .parse::<f64>()
      .map_err(|err| InvalidAmbiguityThresholdText {
        message: format!("'{text}' is not a number: {err}"),
      })?;
    Self::new(value).map_err(|err| InvalidAmbiguityThresholdText {
      message: err.to_string(),
    })
  }
}

impl JsonSchema for AmbiguityThreshold {
  fn inline_schema() -> bool {
    true
  }

  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("AmbiguityThreshold")
  }

  fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
    json_schema!({
      "type": "number",
      "minimum": 0.0,
      "maximum": 1.0,
    })
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Display, Error)]
#[display("the ambiguity threshold must be a probability difference from 0 to 1, got {value}")]
pub struct InvalidAmbiguityThreshold {
  #[error(not(source))]
  pub value: f64,
}

#[derive(Clone, Debug, PartialEq, Eq, Display, Error)]
#[display("{message}")]
pub struct InvalidAmbiguityThresholdText {
  #[error(not(source))]
  message: String,
}
