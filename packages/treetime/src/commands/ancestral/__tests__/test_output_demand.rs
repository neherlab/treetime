#[cfg(test)]
mod tests {
  use crate::commands::ancestral::result::AncestralOutputMaps;
  use crate::commands::ancestral::run::{collect_ancestral_tree_maps, tree_outputs_need_sequences};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::graph::TreeWriteKind;
  use treetime_io::nwk::NwkStyle;

  fn outputs(kinds: &[TreeWriteKind]) -> BTreeMap<TreeWriteKind, PathBuf> {
    kinds.iter().map(|kind| (kind.clone(), PathBuf::from("out"))).collect()
  }

  /// The gather closure returns one non-empty edge-mutation entry so a collected result is
  /// distinguishable from the empty default returned on a skip.
  fn marked_gather(calls: &mut usize) -> Result<AncestralOutputMaps, Report> {
    *calls += 1;
    Ok(AncestralOutputMaps {
      edge_mutations: BTreeMap::from([(GraphEdgeKey(0), vec![])]),
      ..AncestralOutputMaps::default()
    })
  }

  #[test]
  fn test_output_demand_collects_once_for_sequence_writer() -> Result<(), Report> {
    let tree_outputs = outputs(&[TreeWriteKind::nwk(NwkStyle::Plain)]);
    let mut calls = 0;
    let maps = collect_ancestral_tree_maps(&tree_outputs, || marked_gather(&mut calls))?;
    assert_eq!(1, calls, "a sequence-reading writer collects the maps exactly once");
    assert_eq!(
      1,
      maps.edge_mutations.len(),
      "the collected maps are returned to the writer"
    );
    Ok(())
  }

  #[test]
  fn test_output_demand_shares_single_collection_across_writers() -> Result<(), Report> {
    let tree_outputs = outputs(&[
      TreeWriteKind::nwk(NwkStyle::Plain),
      TreeWriteKind::Auspice,
      TreeWriteKind::GraphJson,
    ]);
    let mut calls = 0;
    collect_ancestral_tree_maps(&tree_outputs, || marked_gather(&mut calls))?;
    assert_eq!(1, calls, "several sequence writers share one collection");
    Ok(())
  }

  #[test]
  fn test_output_demand_skips_collection_for_topology_only_writers() -> Result<(), Report> {
    let tree_outputs = outputs(&[TreeWriteKind::GraphJson, TreeWriteKind::Dot]);
    let mut calls = 0;
    let maps = collect_ancestral_tree_maps(&tree_outputs, || marked_gather(&mut calls))?;
    assert_eq!(0, calls, "topology-only writers never expand the sequences");
    assert!(maps.root_sequence.is_none());
    assert!(maps.edge_mutations.is_empty());
    Ok(())
  }

  #[test]
  fn test_output_demand_skips_collection_when_no_tree_writer_selected() -> Result<(), Report> {
    let tree_outputs = outputs(&[]);
    let mut calls = 0;
    collect_ancestral_tree_maps(&tree_outputs, || marked_gather(&mut calls))?;
    assert_eq!(0, calls, "a GTR-only or FASTA-only request skips sequence collection");
    Ok(())
  }

  #[test]
  fn test_output_demand_predicate_matches_writer_needs() {
    assert!(tree_outputs_need_sequences(&outputs(&[TreeWriteKind::nwk(
      NwkStyle::Plain
    )])));
    assert!(tree_outputs_need_sequences(&outputs(&[TreeWriteKind::nexus(
      NwkStyle::Plain
    )])));
    assert!(tree_outputs_need_sequences(&outputs(&[TreeWriteKind::Auspice])));
    assert!(tree_outputs_need_sequences(&outputs(&[TreeWriteKind::MatPb])));
    assert!(tree_outputs_need_sequences(&outputs(&[TreeWriteKind::MatJson])));

    assert!(!tree_outputs_need_sequences(&outputs(&[TreeWriteKind::GraphJson])));
    assert!(!tree_outputs_need_sequences(&outputs(&[TreeWriteKind::Dot])));
    assert!(!tree_outputs_need_sequences(&outputs(&[
      TreeWriteKind::GraphJson,
      TreeWriteKind::Dot,
    ])));
    assert!(!tree_outputs_need_sequences(&outputs(&[])));

    // A topology-only selection mixed with a sequence writer still needs the collection.
    assert!(tree_outputs_need_sequences(&outputs(&[
      TreeWriteKind::GraphJson,
      TreeWriteKind::Auspice,
    ])));
  }
}
