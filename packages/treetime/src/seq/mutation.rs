use crate::alphabet::alphabet::Alphabet;
use crate::seq::indel::{InDel, InDelKind};
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
use treetime_primitives::Seq;
use treetime_utils::error::to_eyre_error;

static NUC_MUT_RE: LazyLock<Regex> = LazyLock::new(|| {
  Regex::new(r"^(?P<ref>[A-Z])(?P<pos>\d{1,10})(?P<qry>[A-Z])$").expect("NUC_MUT_RE regex compilation")
});

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct Mutation {
  pub track: MutationTrack,
  pub event: MutationEvent,
}

impl Mutation {
  pub fn substitution(track: MutationTrack, substitution: Sub) -> Self {
    Self {
      track,
      event: MutationEvent::Substitution(substitution),
    }
  }

  pub fn indel(track: MutationTrack, indel: &InDel) -> Result<Self, Report> {
    let segment = AlignedMutation::new(indel.range, indel.seq.clone())?;
    let event = match indel.kind {
      InDelKind::Insertion => MutationEvent::Insertion(segment),
      InDelKind::Deletion => MutationEvent::Deletion(segment),
    };
    Ok(Self { track, event })
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub enum MutationTrack {
  Nucleotide,
  AminoAcid(String),
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub enum MutationEvent {
  Substitution(Sub),
  Insertion(AlignedMutation),
  Deletion(AlignedMutation),
}

pub fn mutation_event_strings(event: &MutationEvent) -> Result<Vec<String>, Report> {
  match event {
    MutationEvent::Substitution(substitution) => Ok(vec![substitution.to_string()]),
    MutationEvent::Insertion(segment) => segment
      .sequence
      .iter()
      .enumerate()
      .map(|(offset, state)| mutation_position(segment.range.0, offset).map(|position| format!("-{position}{state}")))
      .collect(),
    MutationEvent::Deletion(segment) => segment
      .sequence
      .iter()
      .enumerate()
      .map(|(offset, state)| mutation_position(segment.range.0, offset).map(|position| format!("{state}{position}-")))
      .collect(),
  }
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct AlignedMutation {
  pub range: (usize, usize),
  pub sequence: Seq,
}

impl AlignedMutation {
  pub fn new(range: (usize, usize), sequence: Seq) -> Result<Self, Report> {
    if range.0 >= range.1 {
      return make_error!(
        "Aligned mutation range must be non-empty and ordered, but found {}..{}",
        range.0,
        range.1
      );
    }
    let range_length = range
      .1
      .checked_sub(range.0)
      .ok_or_else(|| eyre::eyre!("Aligned mutation range underflow for {}..{}", range.0, range.1))?;
    if sequence.len() != range_length {
      return make_error!(
        "Aligned mutation range {}..{} has length {range_length}, but its sequence has length {}",
        range.0,
        range.1,
        sequence.len()
      );
    }
    Ok(Self { range, sequence })
  }

  pub fn one_based_inclusive_range(&self) -> Result<(usize, usize), Report> {
    let start = self
      .range
      .0
      .checked_add(1)
      .ok_or_else(|| eyre::eyre!("Mutation start coordinate overflow at {}", self.range.0))?;
    Ok((start, self.range.1))
  }
}

fn mutation_position(start: usize, offset: usize) -> Result<usize, Report> {
  start
    .checked_add(offset)
    .and_then(|position| position.checked_add(1))
    .ok_or_else(|| eyre::eyre!("Mutation coordinate overflow at start {start} and offset {offset}"))
}

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
  debug_assert!(
    parent_subs.windows(2).all(|w| matches!(w, [a, b] if a.pos() < b.pos())),
    "parent_subs not sorted by unique position"
  );
  debug_assert!(
    child_subs.windows(2).all(|w| matches!(w, [a, b] if a.pos() < b.pos())),
    "child_subs not sorted by unique position"
  );

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
        debug_assert_eq!(
          ps.qry(),
          cs.reff(),
          "Substitution chain broken at position {}: parent produces {} but child expects {}",
          ps.pos(),
          ps.qry(),
          cs.reff()
        );
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
          Sub::new(reff, pos, qry)
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
