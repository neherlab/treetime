use derive_more::{Display, Error};
use deser::adapters::TryFromInto;
use deser::{Deserialize, Serialize};
use regex::{Regex, regex};
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use std::borrow::Cow;
use std::fmt;
use std::iter::StepBy;
use std::ops::Range;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[deser(as = TryFromInto<String>)]
pub struct SiteRange {
  first: usize,
  last: usize,
  step: usize,
}

impl SiteRange {
  pub fn last(self) -> usize {
    self.last
  }

  pub fn columns(self) -> StepBy<Range<usize>> {
    ((self.first - 1)..self.last).step_by(self.step)
  }
}

impl fmt::Display for SiteRange {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    let Self { first, last, step } = *self;
    match (first == last, step) {
      (true, 1) => write!(f, "{first}"),
      (false, 1) => write!(f, "{first}-{last}"),
      _ => write!(f, "{first}-{last}\\{step}"),
    }
  }
}

impl FromStr for SiteRange {
  type Err = InvalidSiteRange;

  fn from_str(text: &str) -> Result<Self, InvalidSiteRange> {
    let invalid = |reason: &str| InvalidSiteRange {
      message: format!("site range '{text}': {reason}"),
    };
    let not_a_range = || {
      invalid(
        "write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b \
         and k",
      )
    };
    let captures = site_range_regex().captures(text).ok_or_else(not_a_range)?;
    let number = |index: usize| {
      captures
        .get(index)
        .map(|number| number.as_str().parse::<usize>().map_err(|_err| not_a_range()))
        .transpose()
    };
    let first = number(1)?.ok_or_else(not_a_range)?;
    let last = number(2)?.unwrap_or(first);
    let step = number(3)?.unwrap_or(1);
    if first == 0 {
      return Err(invalid("sites are numbered from 1"));
    }
    if first > last {
      return Err(invalid(&format!(
        "the first site {first} is after the last site {last}"
      )));
    }
    if step == 0 {
      return Err(invalid("the step must be at least 1"));
    }
    Ok(Self { first, last, step })
  }
}

impl TryFrom<String> for SiteRange {
  type Error = InvalidSiteRange;

  fn try_from(text: String) -> Result<Self, InvalidSiteRange> {
    text.parse()
  }
}

impl From<SiteRange> for String {
  fn from(range: SiteRange) -> Self {
    range.to_string()
  }
}

impl JsonSchema for SiteRange {
  fn inline_schema() -> bool {
    true
  }

  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("SiteRange")
  }

  fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
    json_schema!({
      "type": "string",
      "pattern": site_range_regex().as_str(),
    })
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Display, Error)]
#[display("{message}")]
pub struct InvalidSiteRange {
  #[error(not(source))]
  message: String,
}

fn site_range_regex() -> &'static Regex {
  regex!(r"^([0-9]+)(?:-([0-9]+)(?:\\([0-9]+))?)?$")
}
