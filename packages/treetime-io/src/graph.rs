use crate::auspice::auspice_write_file;
use crate::graphviz::{EdgeToGraphviz, NodeToGraphviz, graphviz_write_file};
use crate::nex::{NexWriteOptions, nex_write_file_with};
use crate::nwk::{CommentProviders, EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_file_with};
use crate::phyloxml::{PhyloxmlJsonOptions, phyloxml_json_write_file, phyloxml_write_file};
use crate::tree_ir::auspice::TreeIrAuspiceWriter;
use crate::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode};
use crate::tree_ir::usher::TreeIrUsherWriter;
use crate::usher_mat::{UsherMatJsonOptions, usher_mat_json_write_file, usher_mat_pb_write_file};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::make_report;
use util_newick::NwkStyle;

/// Format-neutral intermediate representation graph for the formats that cannot be
/// derived from the domain graph alone (PhyloXML, Auspice v2, UShER MAT).
pub type TreeIrGraph = Graph<TreeIrNode, TreeIrEdge, TreeIrData>;

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

/// Write the requested tree outputs.
///
/// Newick, Nexus, Graphviz, and graph-JSON are serialized directly from the domain
/// `graph`. PhyloXML, Auspice v2, and UShER MAT are serialized from the
/// format-neutral `ir` graph, which the command builds by projecting its analysis
/// result and partition data. When an IR format is requested but no `ir` graph is
/// provided, the command does not support that format and a clear error is returned.
pub fn write_tree_outputs<N, E, D>(
  graph: &Graph<N, E, D>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
  ir: Option<&TreeIrGraph>,
) -> Result<(), Report>
where
  N: GraphNode + Named + NodeToNwk + NodeToGraphviz + Serialize,
  E: GraphEdge + HasBranchLength + EdgeToNwk + EdgeToGraphviz + Serialize,
  D: Send + Sync + Default + Serialize,
{
  for (kind, path) in outputs {
    match kind {
      TreeWriteKind::Nwk(spec) => {
        let options = NwkWriteOptions {
          style: spec.style,
          ..NwkWriteOptions::default()
        };
        nwk_write_file_with(path, graph, &options, providers)?;
      },
      TreeWriteKind::Nexus(spec) => {
        let options = NexWriteOptions {
          style: spec.style,
          ..NexWriteOptions::default()
        };
        nex_write_file_with(path, graph, &options, providers)?;
      },
      TreeWriteKind::GraphJson => {
        json_write_file(path, graph, JsonPretty(true))?;
      },
      TreeWriteKind::Dot => {
        graphviz_write_file(path, graph)?;
      },
      TreeWriteKind::Auspice => {
        let ir = require_ir(ir, "Auspice v2 JSON")?;
        auspice_write_file::<TreeIrAuspiceWriter, _, _, _>(path, ir)?;
      },
      TreeWriteKind::Phyloxml => {
        let ir = require_ir(ir, "PhyloXML")?;
        phyloxml_write_file(path, ir)?;
      },
      TreeWriteKind::PhyloxmlJson => {
        let ir = require_ir(ir, "PhyloXML JSON")?;
        phyloxml_json_write_file(path, ir, &PhyloxmlJsonOptions::default())?;
      },
      TreeWriteKind::MatPb => {
        let ir = require_ir(ir, "UShER MAT protobuf")?;
        usher_mat_pb_write_file::<TreeIrUsherWriter, _, _, _>(path, ir)?;
      },
      TreeWriteKind::MatJson => {
        let ir = require_ir(ir, "UShER MAT JSON")?;
        usher_mat_json_write_file::<TreeIrUsherWriter, _, _, _>(path, ir, &UsherMatJsonOptions::default())?;
      },
    }
  }
  Ok(())
}

fn require_ir<'a>(ir: Option<&'a TreeIrGraph>, format: &str) -> Result<&'a TreeIrGraph, Report> {
  ir.ok_or_else(|| make_report!("{format} output was requested but is not available for this command"))
}
