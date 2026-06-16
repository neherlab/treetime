//! UShER MAT (mutation-annotated tree) adapter for the format-neutral TreeIR graph.
//!
//! UShER MAT carries only nucleotide substitutions (encoded as integers
//! `0:A, 1:C, 2:G, 3:T`) plus an embedded Newick string for topology, names, and
//! branch lengths. Amino-acid substitutions, indels, dates, and discrete traits
//! have no representation in the format and are dropped with a warning. Ambiguous
//! or non-ACGT nucleotides cannot be encoded and are skipped with a warning.

use crate::tree_ir::mutation::{NUC_GENE, TreeIrSub};
use crate::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode};
use crate::usher_mat::{
  UsherGraphContext, UsherMetadata, UsherMutation, UsherMutationList, UsherRead, UsherTree, UsherTreeContext,
  UsherTreeNode, UsherWrite,
};
use eyre::Report;
use log::warn;
use treetime_graph::graph::Graph;
use treetime_primitives::AsciiChar;
use treetime_utils::make_error;

/// Converter writing a TreeIR graph to a UShER MAT.
pub struct TreeIrUsherWriter {
  /// Reference (root) nucleotide sequence, used to fill `ref_nuc`. Falls back to the
  /// parent nucleotide when the reference is unavailable at a position.
  reference: Vec<AsciiChar>,
  warned_non_nuc: bool,
  warned_indel: bool,
  warned_ambiguous: bool,
}

impl UsherWrite<TreeIrNode, TreeIrEdge, TreeIrData> for TreeIrUsherWriter {
  fn new(graph: &Graph<TreeIrNode, TreeIrEdge, TreeIrData>) -> Result<Self, Report> {
    let reference = graph
      .data()
      .read_arc()
      .root_sequence
      .get(NUC_GENE)
      .map(|seq| seq.bytes().filter_map(|b| AsciiChar::try_new(b).ok()).collect())
      .unwrap_or_default();
    Ok(Self {
      reference,
      warned_non_nuc: false,
      warned_indel: false,
      warned_ambiguous: false,
    })
  }

  fn usher_node_from_graph_components(
    &mut self,
    context: &UsherGraphContext<TreeIrNode, TreeIrEdge, TreeIrData>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report> {
    let node_name = context.node.name.clone().unwrap_or_default();

    let mut mutation = vec![];
    if let Some(edge) = context.edge {
      if !edge.indels.is_empty() && !self.warned_indel {
        warn!("UShER MAT format cannot represent indels; dropping them from output");
        self.warned_indel = true;
      }
      for sub in &edge.mutations {
        if !sub.is_nuc() {
          if !self.warned_non_nuc {
            warn!("UShER MAT format carries only nucleotide substitutions; dropping amino-acid mutations");
            self.warned_non_nuc = true;
          }
          continue;
        }
        let (Some(par_nuc), Some(mut_nuc)) = (nuc_to_int(sub.parent), nuc_to_int(sub.child)) else {
          if !self.warned_ambiguous {
            warn!("UShER MAT format cannot encode non-ACGT nucleotides; skipping ambiguous substitutions");
            self.warned_ambiguous = true;
          }
          continue;
        };
        let ref_nuc = self
          .reference
          .get(sub.position.saturating_sub(1))
          .and_then(|c| nuc_to_int(*c))
          .unwrap_or(par_nuc);
        mutation.push(UsherMutation {
          position: sub.position as i32,
          ref_nuc,
          par_nuc,
          mut_nuc: vec![mut_nuc],
          chromosome: String::new(),
        });
      }
    }

    let node = UsherTreeNode {
      node_name,
      condensed_leaves: vec![],
    };
    let metadata = UsherMetadata {
      clade_annotations: vec![],
    };
    Ok((node, UsherMutationList { mutation }, metadata))
  }
}

/// Converter reading a TreeIR graph from a UShER MAT.
pub struct TreeIrUsherReader;

impl UsherRead<TreeIrNode, TreeIrEdge, TreeIrData> for TreeIrUsherReader {
  fn new(_tree: &UsherTree) -> Result<Self, Report> {
    Ok(Self)
  }

  fn usher_data_to_graph_data(&mut self, tree: &UsherTree) -> Result<TreeIrData, Report> {
    let has_mutations = tree.node_mutations.iter().any(|ml| !ml.mutation.is_empty());
    Ok(TreeIrData {
      has_mutations,
      ..TreeIrData::default()
    })
  }

  fn usher_node_to_graph_components(&mut self, context: &UsherTreeContext) -> Result<(TreeIrNode, TreeIrEdge), Report> {
    let node = TreeIrNode {
      name: context.node.name.clone(),
      ..TreeIrNode::default()
    };

    let mut mutations = vec![];
    for m in &context.node.mutations {
      let Some(&mut_nuc) = m.mut_nuc.first() else {
        continue;
      };
      mutations.push(TreeIrSub {
        gene: NUC_GENE.to_owned(),
        position: m.position as usize,
        parent: int_to_nuc(m.par_nuc)?,
        child: int_to_nuc(mut_nuc)?,
      });
    }

    let edge = TreeIrEdge {
      branch_length: (context.node.branch_length != 0.0).then_some(context.node.branch_length),
      mutations,
      ..TreeIrEdge::default()
    };

    Ok((node, edge))
  }
}

/// Encode a nucleotide as the UShER integer code (`0:A, 1:C, 2:G, 3:T`).
///
/// Returns `None` for ambiguous or non-ACGT characters, which UShER cannot represent.
fn nuc_to_int(c: AsciiChar) -> Option<i32> {
  match char::from(c).to_ascii_uppercase() {
    'A' => Some(0),
    'C' => Some(1),
    'G' => Some(2),
    'T' => Some(3),
    _ => None,
  }
}

/// Decode a UShER integer nucleotide code (`0:A, 1:C, 2:G, 3:T`) to a character.
fn int_to_nuc(code: i32) -> Result<AsciiChar, Report> {
  let ch = match code {
    0 => 'A',
    1 => 'C',
    2 => 'G',
    3 => 'T',
    other => return make_error!("Invalid UShER nucleotide code {other}: expected 0..=3"),
  };
  AsciiChar::try_from_char(ch)
}
