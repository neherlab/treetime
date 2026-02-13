use crate::convert::convert::{ConverterEdge, ConverterNode};
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{NwkWriteOptions, format_weight};

impl NodeToGraphviz for ConverterNode {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

impl EdgeToGraphviz for ConverterEdge {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .branch_length
      .map(|branch_length| format_weight(branch_length, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.branch_length
  }
}
