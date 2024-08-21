use crate::make_error;
use crate::utils::error::to_eyre_error;
use eyre::{Report, WrapErr};
use lazy_static::lazy_static;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct Mut {
  pub pos: usize,
  pub qry: char,
  #[serde(rename = "ref")]
  pub reff: char,
}

impl FromStr for Mut {
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
          let reff = reff.as_str().chars().next().unwrap();
          let pos = parse_pos(pos.as_str()).wrap_err_with(|| format!("When parsing mutation position in '{s}'"))?;
          let qry = qry.as_str().chars().next().unwrap();
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

impl fmt::Display for Mut {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    write!(f, "{}{}{}", self.reff, self.pos + 1, self.qry)
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct InDel {
  pub range: (usize, usize),
  pub seq: Vec<char>,
  pub deletion: bool, // deletion if True, insertion if False
}

impl fmt::Display for InDel {
  fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
    let delta_str = if self.deletion {
      format!("{} -> {}", String::from_iter(&self.seq), "-".repeat(self.seq.len()))
    } else {
      format!("{} -> {}", "-".repeat(self.seq.len()), String::from_iter(&self.seq))
    };
    write!(f, "{}--{}: {}", self.range.0, self.range.1, delta_str)
  }
}
