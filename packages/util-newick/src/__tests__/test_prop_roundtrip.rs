#[cfg(test)]
mod tests {
  use crate::parse::newick_from_string;
  use crate::types::{NewickEdgeData, NewickGraph, NewickNodeData, NewickValue, NewickWriteOptions, NwkStyle};
  use crate::write::newick_to_string;
  use proptest::prelude::*;
  use std::collections::BTreeMap;

  fn arb_branch_length() -> impl Strategy<Value = Option<f64>> {
    prop_oneof![Just(None), (0.0001_f64..100.0).prop_map(Some),]
  }

  fn arb_value() -> impl Strategy<Value = NewickValue> {
    prop_oneof![
      any::<bool>().prop_map(NewickValue::Boolean),
      (0.001_f64..1000.0).prop_map(NewickValue::Number),
      "[A-Za-z]{1,6}".prop_map(NewickValue::String),
    ]
  }

  fn arb_attrs() -> impl Strategy<Value = BTreeMap<String, NewickValue>> {
    proptest::collection::btree_map("[a-z]{1,4}", arb_value(), 0..3)
  }

  fn arb_leaf() -> impl Strategy<Value = (NewickNodeData, NewickEdgeData)> {
    // Branch attrs require a branch length for unambiguous round-trip
    // (without `:`, attrs written after the subtree are parsed as node attrs)
    (arb_branch_length(), arb_attrs()).prop_flat_map(|(bl, branch_attrs)| {
      let branch_attrs = if bl.is_some() { branch_attrs } else { BTreeMap::new() };
      let name_strat = "[A-Za-z][A-Za-z0-9_]{0,8}".prop_map(Some);
      (name_strat, Just(bl), Just(branch_attrs), arb_attrs()).prop_map(|(name, bl, branch_attrs, node_attrs)| {
        (
          NewickNodeData {
            name,
            node_attrs,
            raw_comments: Vec::new(),
            hybrid: None,
            children: Vec::new(),
          },
          NewickEdgeData {
            branch_length: bl,
            branch_attrs,
            raw_comments: Vec::new(),
            is_acceptor: false,
          },
        )
      })
    })
  }

  fn arb_tree() -> impl Strategy<Value = NewickGraph> {
    proptest::collection::vec(arb_leaf(), 2..7).prop_map(|leaves| {
      let mut g = NewickGraph::new();
      let root = g.add_node(NewickNodeData::new());
      g.root = root;
      for (i, (mut node_data, edge_data)) in leaves.into_iter().enumerate() {
        // Ensure unique names for deterministic order-insensitive comparison
        node_data.name = Some(format!("T{i}_{}", node_data.name.unwrap_or_default()));
        let child = g.add_node(node_data);
        g.add_edge(root, child, edge_data);
      }
      g
    })
  }

  proptest! {
    #[test]
    fn test_prop_roundtrip_beast(g in arb_tree()) {
      let opts = NewickWriteOptions {
        style: NwkStyle::Beast,
        significant_digits: None,
        decimal_digits: None,
      };
      let written = newick_to_string(&g, &opts).unwrap();
      let parsed = newick_from_string(&written).unwrap();
      prop_assert_eq!(&g, &parsed, "Round-trip failed.\nWritten: {}", written);
    }

    #[test]
    fn test_prop_roundtrip_nhx(g in arb_tree()) {
      // NHX converts typed values to strings, so we verify topology and branch lengths only
      let opts = NewickWriteOptions {
        style: NwkStyle::Nhx,
        significant_digits: None,
        decimal_digits: None,
      };
      let written = newick_to_string(&g, &opts).unwrap();
      let parsed = newick_from_string(&written).unwrap();
      // NHX round-trip: topology, names, branch lengths match.
      // Attrs may differ in type (Number -> String) so we only check topology.
      prop_assert_eq!(g.nodes.len(), parsed.nodes.len());
      prop_assert_eq!(g.edges.len(), parsed.edges.len());
      for (e1, e2) in g.edges.iter().zip(parsed.edges.iter()) {
        match (e1.data.branch_length, e2.data.branch_length) {
          (Some(a), Some(b)) => {
            prop_assert!((a - b).abs() < 1e-10,
              "Branch length mismatch: {a} vs {b}");
          }
          (None, None) => {}
          _ => prop_assert!(false, "Branch length presence mismatch"),
        }
      }
    }

    #[test]
    fn test_prop_idempotent_write_beast(g in arb_tree()) {
      let opts = NewickWriteOptions {
        style: NwkStyle::Beast,
        significant_digits: None,
        decimal_digits: None,
      };
      let w1 = newick_to_string(&g, &opts).unwrap();
      let parsed = newick_from_string(&w1).unwrap();
      let w2 = newick_to_string(&parsed, &opts).unwrap();
      prop_assert_eq!(w1, w2, "Write not idempotent");
    }

    #[test]
    fn test_prop_idempotent_write_plain(g in arb_tree()) {
      let opts = NewickWriteOptions {
        style: NwkStyle::Plain,
        significant_digits: None,
        decimal_digits: None,
      };
      let w1 = newick_to_string(&g, &opts).unwrap();
      let parsed = newick_from_string(&w1).unwrap();
      let w2 = newick_to_string(&parsed, &opts).unwrap();
      prop_assert_eq!(w1, w2, "Write not idempotent");
    }
  }
}
