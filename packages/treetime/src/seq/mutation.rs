use crate::make_error;
use crate::representation::seq_char::AsciiChar;
use crate::utils::error::to_eyre_error;
use eyre::{Report, WrapErr};
use lazy_static::lazy_static;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct Sub {
  pub pos: usize,
  pub qry: AsciiChar,
  #[serde(rename = "ref")]
  pub reff: AsciiChar,
}

impl FromStr for Sub {
  type Err = Report;

  /// Parses nucleotide substitution from string. Expects IUPAC notation commonly used in bioinformatics.
  fn from_str(s: &str) -> Result<Self, Self::Err> {
    const NUC_MUT_REGEX: &str = r"((?P<ref>[A-Z-])(?P<pos>\d{1,10})(?P<qry>[A-Z-]))";
    lazy_static! {
      static ref RE: Regex = Regex::new(NUC_MUT_REGEX)
        .wrap_err_with(|| format!("When compiling regular expression '{NUC_MUT_REGEX}'"))
        .unwrap();
    }
    if let Some(captures) = RE.captures(s) {
      return match (captures.name("ref"), captures.name("pos"), captures.name("qry")) {
        (Some(reff), Some(pos), Some(qry)) => {
          let reff = AsciiChar(reff.as_str().bytes().next().unwrap());
          let pos = parse_pos(pos.as_str()).wrap_err_with(|| format!("When parsing mutation position in '{s}'"))?;
          let qry = AsciiChar(qry.as_str().bytes().next().unwrap());
          Ok(Self { pos, qry, reff })
        }
        _ => make_error!("Unable to parse nucleotide mutation: '{s}'"),
      };
    }
    make_error!("Unable to parse nucleotide mutation: '{s}'")
  }
}

/// Parses position from string.
///
/// Note that Nextclade uses 0-based indices, including for positions in sequences. However, in bioinformatics it is
/// more common to see 1-based indexing. We perform the conversion here.
pub fn parse_pos(s: &str) -> Result<usize, Report> {
  let pos = to_eyre_error(s.parse::<usize>()).wrap_err_with(|| format!("Unable to parse position: '{s}'"))?;
  if pos < 1 {
    return make_error!("Mutation position is expected to be >= 1");
  }
  Ok(pos - 1)
}

impl fmt::Display for Sub {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{}{}{}", self.reff, self.pos + 1, self.qry)
  }
}
