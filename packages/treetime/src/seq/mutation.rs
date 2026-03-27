use crate::alphabet::alphabet::Alphabet;
use crate::{make_error, make_internal_error};
use eyre::{Report, WrapErr};
use getset::CopyGetters;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fmt;
use std::str::FromStr;
use std::sync::LazyLock;
use treetime_primitives::AsciiChar;
use treetime_utils::error::to_eyre_error;

static NUC_MUT_RE: LazyLock<Regex> = LazyLock::new(|| {
  Regex::new(r"((?P<ref>[A-Z-])(?P<pos>\d{1,10})(?P<qry>[A-Z-]))").expect("NUC_MUT_RE regex compilation")
});

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq, CopyGetters)]
#[getset(get_copy = "pub")]
pub struct Sub {
  pos: usize,
  qry: AsciiChar,
  #[serde(rename = "ref")]
  reff: AsciiChar,
}

impl Sub {
  pub fn new<P: Into<usize>>(reff: AsciiChar, pos: P, qry: AsciiChar) -> Result<Self, Report> {
    let pos = pos.into();

    if qry == AsciiChar::from_byte_unchecked(b'-') || reff == AsciiChar::from_byte_unchecked(b'-') {
      return make_internal_error!("Substitution cannot be from or to gap, but found: '{reff}{pos}{qry}'");
    }

    Ok(Self { pos, qry, reff })
  }

  pub fn check_determined(&self, alphabet: &Alphabet) -> Result<(), Report> {
    if !alphabet.is_determined(self.qry()) || !alphabet.is_determined(self.reff()) {
      make_internal_error!("Substitution is not determined: '{self}'")
    } else {
      Ok(())
    }
  }

  pub fn check_canonical(&self, alphabet: &Alphabet) -> Result<(), Report> {
    if !alphabet.is_canonical(self.qry()) || !alphabet.is_canonical(self.reff()) {
      make_internal_error!("Substitution is not canonical: '{self}'")
    } else {
      Ok(())
    }
  }

  /// Invert the substitution direction by swapping ref and query
  pub fn invert(&mut self) {
    std::mem::swap(&mut self.reff, &mut self.qry);
  }
}

/// Compose substitutions from two consecutive edges into net substitutions.
///
/// When collapsing node B between parent A and child C, substitutions on
/// edges A→B (`parent_subs`) and B→C (`child_subs`) are composed to produce
/// the net A→C substitutions:
///
/// - Non-overlapping positions: kept as-is from whichever edge
/// - Same position, chain: parent ref→X + child X→qry = net ref→qry
/// - Same position, cancellation: parent ref→X + child X→ref = no net change
///
/// Both input slices must be sorted by position with at most one entry per
/// position. The output is sorted by position.
pub fn compose_substitutions(parent_subs: &[Sub], child_subs: &[Sub]) -> Result<Vec<Sub>, Report> {
  let mut result = Vec::with_capacity(parent_subs.len() + child_subs.len());
  let mut pi = 0;
  let mut ci = 0;

  while pi < parent_subs.len() && ci < child_subs.len() {
    let ps = &parent_subs[pi];
    let cs = &child_subs[ci];

    match ps.pos().cmp(&cs.pos()) {
      Ordering::Less => {
        result.push(ps.clone());
        pi += 1;
      },
      Ordering::Greater => {
        result.push(cs.clone());
        ci += 1;
      },
      Ordering::Equal => {
        // Compose: net change is parent.reff → child.qry
        if ps.reff() != cs.qry() {
          result.push(Sub::new(ps.reff(), ps.pos(), cs.qry())?);
        }
        // parent.reff == child.qry: mutations cancel, no net change
        pi += 1;
        ci += 1;
      },
    }
  }

  result.extend_from_slice(&parent_subs[pi..]);
  result.extend_from_slice(&child_subs[ci..]);

  Ok(result)
}

impl FromStr for Sub {
  type Err = Report;

  /// Parses nucleotide substitution from string. Expects IUPAC notation commonly used in bioinformatics.
  fn from_str(s: &str) -> Result<Self, Self::Err> {
    if let Some(captures) = NUC_MUT_RE.captures(s) {
      return match (captures.name("ref"), captures.name("pos"), captures.name("qry")) {
        (Some(reff), Some(pos), Some(qry)) => {
          let reff = AsciiChar::try_new(reff.as_str().bytes().next().unwrap())
            .wrap_err_with(|| format!("When parsing ref character in '{s}'"))?;
          let pos = parse_pos(pos.as_str()).wrap_err_with(|| format!("When parsing mutation position in '{s}'"))?;
          let qry = AsciiChar::try_new(qry.as_str().bytes().next().unwrap())
            .wrap_err_with(|| format!("When parsing qry character in '{s}'"))?;
          Ok(Self { pos, qry, reff })
        },
        _ => make_error!("Unable to parse nucleotide mutation: '{s}'"),
      };
    }
    make_error!("Unable to parse nucleotide mutation: '{s}'")
  }
}

/// Parse position from 1-based bioinformatics notation to 0-based index.
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
