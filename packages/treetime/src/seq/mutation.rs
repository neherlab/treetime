use crate::alphabet::alphabet::Alphabet;
use crate::seq::indel::{InDel, InDelKind};
use crate::{make_error, make_internal_error};
use derive_more::Display;
use eyre::{Report, WrapErr};
use getset::CopyGetters;
use regex::regex;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::str::FromStr;
use treetime_primitives::AsciiChar;
use treetime_primitives::Seq;
use treetime_utils::error::to_eyre_error;

pub(crate) fn combine_edge_mutations(
  subs: Vec<Sub>,
  indels: &[InDel],
  track: &MutationTrack,
) -> Result<Vec<Mutation>, Report> {
  subs
    .into_iter()
    .map(|substitution| Ok(Mutation::substitution(track.clone(), substitution)))
    .chain(indels.iter().map(|indel| Mutation::indel(track.clone(), indel)))
    .collect()
}

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
#[serde(rename_all = "kebab-case")]
pub enum MutationTrack {
  Nucleotide,
  AminoAcid(String),
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
#[serde(rename_all = "kebab-case")]
pub enum MutationEvent {
  Substitution(Sub),
  Insertion(AlignedMutation),
  Deletion(AlignedMutation),
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct AlignedMutation {
  pub range: (usize, usize),
  pub sequence: Seq,
}

impl AlignedMutation {
  pub(crate) fn new(range: (usize, usize), sequence: Seq) -> Result<Self, Report> {
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
}

fn mutation_position(start: usize, offset: usize) -> Result<usize, Report> {
  start
    .checked_add(offset)
    .and_then(|position| position.checked_add(1))
    .ok_or_else(|| eyre::eyre!("Mutation coordinate overflow at start {start} and offset {offset}"))
}

pub(crate) fn compose_substitutions(parent_subs: &[Sub], child_subs: &[Sub]) -> Result<Vec<Sub>, Report> {
  debug_assert!(
    parent_subs.is_sorted_by(|a, b| a.pos() < b.pos()),
    "parent_subs not sorted by unique position"
  );
  debug_assert!(
    child_subs.is_sorted_by(|a, b| a.pos() < b.pos()),
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
        if ps.reff() != cs.qry() {
          result.push(Sub::new(ps.reff(), ps.pos(), cs.qry())?);
        }
        pi += 1;
        ci += 1;
      },
    }
  }

  result.extend_from_slice(&parent_subs[pi..]);
  result.extend_from_slice(&child_subs[ci..]);

  Ok(result)
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq, CopyGetters, Display)]
#[getset(get_copy = "pub")]
#[display("{reff}{}{qry}", pos + 1)]
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

  pub(crate) fn check_determined(&self, alphabet: &Alphabet) -> Result<(), Report> {
    if !alphabet.is_determined(self.qry()) || !alphabet.is_determined(self.reff()) {
      make_internal_error!("Substitution is not determined: '{self}'")
    } else {
      Ok(())
    }
  }

  pub(crate) fn check_canonical(&self, alphabet: &Alphabet) -> Result<(), Report> {
    if !alphabet.is_canonical(self.qry()) || !alphabet.is_canonical(self.reff()) {
      make_internal_error!("Substitution is not canonical: '{self}'")
    } else {
      Ok(())
    }
  }

  pub(crate) fn invert(&mut self) {
    std::mem::swap(&mut self.reff, &mut self.qry);
  }
}

impl FromStr for Sub {
  type Err = Report;

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  fn from_str(s: &str) -> Result<Self, Self::Err> {
    if let Some(captures) = regex!(r"^(?P<ref>[A-Z])(?P<pos>\d{1,10})(?P<qry>[A-Z])$").captures(s) {
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

fn parse_pos(s: &str) -> Result<usize, Report> {
  let pos = to_eyre_error(s.parse::<usize>()).wrap_err_with(|| format!("Unable to parse position: '{s}'"))?;
  if pos < 1 {
    return make_error!("Mutation position is expected to be >= 1");
  }
  Ok(pos - 1)
}
