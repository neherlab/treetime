use crate::graphviz::{EdgeToGraphviz, NodeToGraphviz, graphviz_write_file};
use crate::nex::{NexWriteOptions, nex_write_file_with};
use crate::nwk::{CommentProviders, EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_file_with};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::{make_error, make_report};
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

pub fn write_tree_outputs<N, E, D>(
  graph: &Graph<N, E, D>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
  auspice_writer: Option<&dyn AuspiceWriter<N, E, D>>,
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
        let writer = auspice_writer.ok_or_else(|| {
          make_report!(
            "Auspice output requested but no auspice writer provided. \
             Auspice JSON requires domain-specific data only available on certain commands."
          )
        })?;
        writer.write_auspice(graph, path)?;
      },
      TreeWriteKind::Phyloxml | TreeWriteKind::PhyloxmlJson => {
        return make_error!("PhyloXML output is not yet implemented for analysis commands");
      },
      TreeWriteKind::MatPb | TreeWriteKind::MatJson => {
        return make_error!("UShER MAT output is not yet implemented for analysis commands");
      },
    }
  }
  Ok(())
}

/// Trait for command-specific auspice JSON writing.
///
/// Different commands have different data to include in the auspice JSON
/// (confidence intervals, mutation counts, etc). This trait abstracts
/// the write operation so `write_tree_outputs` can dispatch without
/// depending on command-specific types.
pub trait AuspiceWriter<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  fn write_auspice(&self, graph: &Graph<N, E, D>, path: &Path) -> Result<(), Report>;
}
