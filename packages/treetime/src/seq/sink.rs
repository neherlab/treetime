use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

/// Streaming sink for reconstructed per-node sequences, keyed by graph node.
///
/// The reconstruction drivers emit each node's sequence by graph key as it is produced, so the whole
/// set never resides in memory at once. Core emits keys only; a caller-supplied sink maps each key to an
/// output name and description and encodes the record. [`SeqSink::on_topology`] is called once with the
/// final tree before the first [`SeqSink::emit`], so a sink can resolve names against the topology a
/// late reroot or polytomy resolution produced: the synthetic names core assigns are not written onto
/// the graph, so a sink reproduces them from the final topology rather than reading them off the nodes.
pub trait SeqSink {
  /// Announce the final tree before any sequence is emitted.
  fn on_topology(&mut self, graph: &Graph) -> Result<(), Report>;
  /// Emit one node's reconstructed sequence for one track.
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report>;
}

/// One reconstructed sequence handed to a [`SeqSink`]: the graph node it belongs to, the track it was
/// reconstructed on, and the sequence itself.
pub struct SeqItem<'a> {
  pub key: GraphNodeKey,
  pub track: SeqTrack<'a>,
  pub seq: &'a Seq,
}

/// The track a reconstructed sequence belongs to: the nucleotide partition, or a named amino-acid CDS.
pub enum SeqTrack<'a> {
  Nuc,
  Aa(&'a str),
}
