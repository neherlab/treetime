use util_newick::NwkStyle;

/// Dispatch tag for tree output formats. Used as keys in the resolved output map.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum TreeWriteKind {
  Nwk(NwkWriteSpec),
  Nexus(NwkWriteSpec),
  Auspice,
  Phyloxml,
  PhyloxmlJson,
  MatPb,
  MatJson,
  GraphJson,
  Dot,
}

impl TreeWriteKind {
  pub fn nwk(style: NwkStyle) -> Self {
    Self::Nwk(NwkWriteSpec { style })
  }

  pub fn nexus(style: NwkStyle) -> Self {
    Self::Nexus(NwkWriteSpec { style })
  }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct NwkWriteSpec {
  pub style: NwkStyle,
}
