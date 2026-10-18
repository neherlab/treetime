use crate::seq::mutation::MutationTrack;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

pub trait SeqSink {
  fn on_topology(&mut self, _graph: &Graph) -> Result<(), Report> {
    Ok(())
  }

  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report>;
}

pub struct SeqItem<'a> {
  pub key: GraphNodeKey,
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
