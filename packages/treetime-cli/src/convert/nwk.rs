use crate::convert::convert::{ConverterEdge, ConverterNode};
use eyre::Report;
use std::collections::BTreeMap;
use treetime::io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk};

impl NodeFromNwk for ConverterNode {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|n| n.as_ref().to_owned()),
    })
  }
}

impl NodeToNwk for ConverterNode {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

impl EdgeFromNwk for ConverterEdge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self { weight })
  }
}

impl EdgeToNwk for ConverterEdge {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight
  }
}
