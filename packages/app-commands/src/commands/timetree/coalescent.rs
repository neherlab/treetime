use derive_more::{Display, Error};
use deser::adapters::TryFromInto;
use deser::{Deserialize, Serialize};
use deser_serde::Serde;
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use serde_json::Value;
use std::borrow::Cow;
use std::str::FromStr;
use treetime::timetree::coalescent_timescale::{CoalescentMode, FixedTimescale};

const DAYS_PER_YEAR: f64 = 365.25;

const WEEKS_PER_YEAR: f64 = DAYS_PER_YEAR / 7.0;

const FIXED_TIMESCALE_PATTERN: &str =
  r"^\+?([0-9]*[1-9][0-9]*\.?[0-9]*|[0-9]*\.[0-9]*[1-9][0-9]*)([eE][+-]?[0-9]+)?[ywd]?$";

pub fn parse_coalescent(text: &str) -> Result<CoalescentArg, InvalidCoalescent> {
  match text {
    "opt" => Ok(CoalescentArg::Opt),
    "skyline" => Ok(CoalescentArg::Skyline),
    _ => parse_fixed_timescale(text)
      .map(CoalescentArg::Fixed)
      .ok_or_else(|| InvalidCoalescent { text: text.to_owned() }),
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Display, Serialize, Deserialize)]
#[deser(as = TryFromInto<CoalescentValue>)]
pub enum CoalescentArg {
  #[display("{_0}")]
  Fixed(FixedTimescaleArg),
  #[display("opt")]
  Opt,
  #[display("skyline")]
  Skyline,
}

impl CoalescentArg {
  pub fn mode(self) -> CoalescentMode {
    match self {
      Self::Fixed(fixed) => CoalescentMode::Fixed(fixed.timescale()),
      Self::Opt => CoalescentMode::Constant,
      Self::Skyline => CoalescentMode::Skyline,
    }
  }
}

impl FromStr for CoalescentArg {
  type Err = InvalidCoalescent;

  fn from_str(text: &str) -> Result<Self, InvalidCoalescent> {
    parse_coalescent(text)
  }
}

impl TryFrom<CoalescentValue> for CoalescentArg {
  type Error = InvalidCoalescent;

  fn try_from(value: CoalescentValue) -> Result<Self, InvalidCoalescent> {
    match value {
      CoalescentValue::Number(number) => FixedTimescaleArg::new(number, TimescaleUnit::SubstitutionsPerSite)
        .map(Self::Fixed)
        .ok_or_else(|| InvalidCoalescent {
          text: number.to_string(),
        }),
      CoalescentValue::Text(text) => parse_coalescent(&text),
      CoalescentValue::Other(value) => Err(InvalidCoalescent {
        text: value.to_string(),
      }),
    }
  }
}

impl From<CoalescentArg> for CoalescentValue {
  fn from(arg: CoalescentArg) -> Self {
    match arg {
      CoalescentArg::Fixed(FixedTimescaleArg {
        value,
        unit: TimescaleUnit::SubstitutionsPerSite,
      }) => Self::Number(value),
      arg => Self::Text(arg.to_string()),
    }
  }
}

impl JsonSchema for CoalescentArg {
  fn inline_schema() -> bool {
    true
  }

  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("CoalescentArg")
  }

  fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
    json_schema!({
      "oneOf": [
        {
          "type": "number",
          "exclusiveMinimum": 0.0,
          "description": "Fixed time scale in substitutions per site",
        },
        {
          "type": "string",
          "pattern": FIXED_TIMESCALE_PATTERN,
          "description": "Fixed time scale: a number in substitutions per site, or a number directly followed by the time unit y, w or d",
        },
        {
          "type": "string",
          "enum": ["opt", "skyline"],
        },
      ],
    })
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Display)]
#[display("{value}{unit}")]
pub struct FixedTimescaleArg {
  pub value: f64,
  pub unit: TimescaleUnit,
}

impl FixedTimescaleArg {
  fn new(value: f64, unit: TimescaleUnit) -> Option<Self> {
    let is_positive = |number: f64| number.is_finite() && number > 0.0;
    (is_positive(value) && is_positive(unit.in_scale_unit(value))).then_some(Self { value, unit })
  }

  pub fn timescale(self) -> FixedTimescale {
    let scale = self.unit.in_scale_unit(self.value);
    match self.unit {
      TimescaleUnit::SubstitutionsPerSite => FixedTimescale::SubstitutionsPerSite(scale),
      TimescaleUnit::Years | TimescaleUnit::Weeks | TimescaleUnit::Days => FixedTimescale::Years(scale),
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Display)]
pub enum TimescaleUnit {
  #[display("")]
  SubstitutionsPerSite,
  #[display("y")]
  Years,
  #[display("w")]
  Weeks,
  #[display("d")]
  Days,
}

impl TimescaleUnit {
  fn from_suffix(suffix: char) -> Option<Self> {
    match suffix {
      'y' => Some(Self::Years),
      'w' => Some(Self::Weeks),
      'd' => Some(Self::Days),
      _ => None,
    }
  }

  fn in_scale_unit(self, value: f64) -> f64 {
    match self {
      Self::SubstitutionsPerSite | Self::Years => value,
      Self::Weeks => value / WEEKS_PER_YEAR,
      Self::Days => value / DAYS_PER_YEAR,
    }
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Display, Error)]
#[display(
  "'{text}' is not a coalescent prior. Give a fixed time scale as a positive number in substitutions per site (0.01, \
   1e-2) or as a positive number directly followed by the time unit y, w or d (10y, 26w, 90d); opt to fit a constant \
   time scale; or skyline to fit a piecewise-constant one"
)]
pub struct InvalidCoalescent {
  #[error(not(source))]
  pub text: String,
}

fn parse_fixed_timescale(text: &str) -> Option<FixedTimescaleArg> {
  let suffix = text.chars().last()?;
  let (number, unit) = TimescaleUnit::from_suffix(suffix).map_or((text, TimescaleUnit::SubstitutionsPerSite), |unit| {
    (text.strip_suffix(suffix).unwrap_or(text), unit)
  });
  FixedTimescaleArg::new(number.parse().ok()?, unit)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[deser(untagged)]
enum CoalescentValue {
  Number(f64),
  Text(String),
  Other(#[deser(as = Serde)] Value),
}
