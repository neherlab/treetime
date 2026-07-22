use crate::auspice::{AuspiceWrite, auspice_write_file};
use crate::graphviz::{EdgeToGraphviz, NodeToGraphviz, graphviz_write_file};
use crate::nex::{NexWriteOptions, nex_write_file_with};
use crate::nwk::{CommentProviders, EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_file_with};
use crate::phyloxml::{PhyloxmlFromGraph, PhyloxmlJsonOptions, phyloxml_json_write_file, phyloxml_write_file};
use crate::usher_mat::{UsherMatJsonOptions, UsherWrite, usher_mat_json_write_file, usher_mat_pb_write_file};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_utils::io::json::{JsonPretty, json_write_file};
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

/// Write every requested tree format from one authoritative graph.
pub fn write_tree_outputs<C, N, E, D>(
  graph: &Graph<N, E, D>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report>
where
  N: GraphNode + Named + NodeToNwk + NodeToGraphviz + Serialize,
  E: GraphEdge + HasBranchLength + EdgeToNwk + EdgeToGraphviz + Serialize,
  D: Send + Sync + Serialize,
  C: AuspiceWrite<N, E, D> + PhyloxmlFromGraph<N, E, D> + UsherWrite<N, E, D>,
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
        auspice_write_file::<C, _, _, _>(path, graph)?;
      },
      TreeWriteKind::Phyloxml => {
        phyloxml_write_file::<C, _, _, _>(path, graph)?;
      },
      TreeWriteKind::PhyloxmlJson => {
        phyloxml_json_write_file::<C, _, _, _>(path, graph, &PhyloxmlJsonOptions::default())?;
      },
      TreeWriteKind::MatPb => {
        usher_mat_pb_write_file::<C, _, _, _>(path, graph)?;
      },
      TreeWriteKind::MatJson => {
        usher_mat_json_write_file::<C, _, _, _>(path, graph, &UsherMatJsonOptions::default())?;
      },
    }
  }
  Ok(())
}
