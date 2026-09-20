use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

pub trait SeqSink {
  fn on_topology(&mut self, graph: &Graph) -> Result<(), Report>;
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report>;
}

pub struct SeqItem<'a> {
  pub key: GraphNodeKey,
  pub track: SeqTrack<'a>,
  pub seq: &'a Seq,
}

pub enum SeqTrack<'a> {
  Nuc,
  Aa(&'a str),
}
