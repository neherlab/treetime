use crate::convert::convert::{ConverterEdge, ConverterNode};
use treetime::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use treetime::io::nwk::{format_weight, NwkWriteOptions};

impl NodeToGraphviz for ConverterNode {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

impl EdgeToGraphViz for ConverterEdge {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight
  }
}
