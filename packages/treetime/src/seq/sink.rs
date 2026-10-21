use crate::seq::mutation::MutationTrack;
use eyre::Report;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::Seq;

pub trait SeqSink {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report>;
}

pub struct SeqItem<'a> {
  pub key: GraphNodeKey,
  pub name: &'a str,
  pub track: SeqTrack<'a>,
  pub seq: &'a Seq,
  pub emitted: bool,
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

  pub(crate) fn emit(
    &mut self,
    key: GraphNodeKey,
    track: SeqTrack<'_>,
    seq: &Seq,
    emitted: bool,
  ) -> Result<(), Report> {
    self.sink.emit(SeqItem {
      key,
      name: &self.names[key],
      track,
      seq,
      emitted,
    })
  }
}
