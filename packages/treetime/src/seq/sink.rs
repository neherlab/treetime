use crate::partition::marginal::posterior_states::AmbiguityThreshold;
use crate::seq::mutation::MutationTrack;
use eyre::Report;
use std::borrow::Cow;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::{AsciiChar, Seq};

pub trait SeqSink {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report>;

  fn ambiguity(&self) -> Option<AmbiguityThreshold> {
    None
  }
}

pub type AmbiguousSites = BTreeMap<usize, AsciiChar>;

pub struct SeqItem<'a> {
  pub key: GraphNodeKey,
  pub name: &'a str,
  pub track: SeqTrack<'a>,
  pub seq: &'a Seq,
  pub ambiguous_sites: &'a AmbiguousSites,
  pub emitted: bool,
}

impl<'a> SeqItem<'a> {
  pub fn seq_with_ambiguity(&self) -> Cow<'a, Seq> {
    if self.ambiguous_sites.is_empty() {
      return Cow::Borrowed(self.seq);
    }
    let mut seq = self.seq.clone();
    for (&pos, &code) in self.ambiguous_sites {
      seq[pos] = code;
    }
    Cow::Owned(seq)
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqTrack<'a> {
  Nuc(usize),
  Aa(&'a str),
}

impl SeqTrack<'_> {
  pub(crate) fn mutation_track(self) -> MutationTrack {
    match self {
      Self::Nuc(_) => MutationTrack::Nucleotide,
      Self::Aa(cds) => MutationTrack::AminoAcid(cds.to_owned()),
    }
  }
}

pub struct SeqStream<'a> {
  sink: &'a mut dyn SeqSink,
  names: &'a NodeNames,
}

impl<'a> SeqStream<'a> {
  pub fn new(sink: &'a mut dyn SeqSink, names: &'a NodeNames) -> Self {
    Self { sink, names }
  }

  pub(crate) fn reborrow(&mut self) -> SeqStream<'_> {
    SeqStream {
      sink: &mut *self.sink,
      names: self.names,
    }
  }

  pub(crate) fn ambiguity(&self) -> Option<AmbiguityThreshold> {
    self.sink.ambiguity()
  }

  pub(crate) fn emit(
    &mut self,
    key: GraphNodeKey,
    track: SeqTrack<'_>,
    seq: &Seq,
    ambiguous_sites: &AmbiguousSites,
    emitted: bool,
  ) -> Result<(), Report> {
    self.sink.emit(SeqItem {
      key,
      name: &self.names[key],
      track,
      seq,
      ambiguous_sites,
      emitted,
    })
  }
}
