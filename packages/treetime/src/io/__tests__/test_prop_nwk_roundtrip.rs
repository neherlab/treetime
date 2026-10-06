#[cfg(test)]
mod tests {
  use eyre::Report;
  use generators::arb_tree;
  use helpers::{build_graph, read_back};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use std::thread;
  use treetime_graph::tree_view::TreeView;
  use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};

  proptest! {
    #[test]
    fn test_prop_nwk_roundtrip_names_and_lengths(tree in arb_tree()) {
      let (expected, written) = build_graph(&tree).unwrap();

      let actual = read_back(&written).unwrap();

      prop_assert_eq!(expected, actual, "written: {}", written);
    }
  }

  #[test]
  fn test_prop_nwk_roundtrip_deep_tree_on_small_stack() -> Result<(), Report> {
    let depth = 100_000;
    let mut input = "(".repeat(depth);
    input.push_str("A:0.01");
    input.push_str(&",L:0.01):0.01".repeat(depth));
    input.push(';');

    let written = thread::Builder::new()
      .stack_size(2 << 20)
      .spawn(move || -> Result<String, Report> {
        let parse = nwk_read(input.as_bytes())?;
        let names = parse.names();
        nwk_write_str(
          &TreeView::new(&parse.graph)?,
          &names,
          &parse.branch_lengths,
          &NwkWriteOptions::default(),
          &NwkNodeComments::new(),
        )
      })?
      .join()
      .unwrap()?;

    assert_eq!(2 * depth + 1, written.matches(':').count() + 1);
    Ok(())
  }

  mod helpers {
    use super::generators::GenTree;
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::tree_view::TreeView;
    use treetime_io::nwk::{NwkNodeComments, NwkWriteOptions, nwk_read, nwk_write_str};

    pub(super) type NameParentLength = BTreeMap<String, (Option<String>, Option<u64>)>;

    pub(super) fn build_graph(tree: &GenTree) -> Result<(NameParentLength, String), Report> {
      let mut graph = Graph::new();
      let mut names = BTreeMap::new();
      let mut weights = BTreeMap::new();
      let mut expected = BTreeMap::new();
      let mut stack: Vec<(&GenTree, Option<(GraphNodeKey, &str)>, Option<f64>)> = vec![(tree, None, None)];
      while let Some((node, parent, length)) = stack.pop() {
        let key = graph.add_node();
        names.insert(key, Some(node.name.clone()));
        expected.insert(
          node.name.clone(),
          (parent.map(|(_, name)| name.to_owned()), length.map(f64::to_bits)),
        );
        if let Some((parent_key, _)) = parent {
          let edge = graph.add_edge(parent_key, key)?;
          weights.insert(edge, length);
        }
        for (child, child_length) in node.children.iter().rev() {
          stack.push((child, Some((key, &node.name)), *child_length));
        }
      }
      graph.build()?;
      let options = NwkWriteOptions {
        weight_significant_digits: Some(17),
        ..NwkWriteOptions::default()
      };
      let written = nwk_write_str(
        &TreeView::new(&graph)?,
        &names,
        &weights,
        &options,
        &NwkNodeComments::new(),
      )?;
      Ok((expected, written))
    }

    pub(super) fn read_back(written: &str) -> Result<NameParentLength, Report> {
      let parse = nwk_read(written.as_bytes())?;
      let names = parse.names();
      let mut actual = BTreeMap::new();
      for node in parse.graph.get_nodes() {
        let name = names[&node.key()].clone().unwrap_or_default();
        let inbound = node.inbound();
        let (parent, length) = match inbound.first() {
          Some(&edge_key) => {
            let edge = parse.graph.get_edge(edge_key).unwrap();
            (
              names[&edge.source()].clone(),
              parse.branch_lengths[&edge_key].map(f64::to_bits),
            )
          },
          None => (None, None),
        };
        actual.insert(name, (parent, length));
      }
      Ok(actual)
    }
  }

  mod generators {
    use proptest::collection::vec;
    use proptest::prelude::*;

    #[derive(Clone, Debug)]
    pub(super) struct GenTree {
      pub(super) name: String,
      pub(super) children: Vec<(GenTree, Option<f64>)>,
    }

    pub(super) fn arb_tree() -> impl Strategy<Value = GenTree> {
      let leaf = "\\PC{0,6}".prop_map(|name| GenTree {
        name,
        children: Vec::new(),
      });
      let shape = leaf.prop_recursive(5, 48, 4, |inner| {
        ("\\PC{0,6}", vec((inner, proptest::option::of(arb_length())), 1..4))
          .prop_map(|(name, children)| GenTree { name, children })
      });
      shape.prop_map(|mut tree| {
        let mut counter = 0;
        make_names_unique(&mut tree, &mut counter);
        tree
      })
    }

    fn make_names_unique(tree: &mut GenTree, counter: &mut usize) {
      let separator = if tree.children.is_empty() { "" } else { "~" };
      tree.name = format!("{}{separator}{counter}", tree.name);
      *counter += 1;
      for (child, _) in &mut tree.children {
        make_names_unique(child, counter);
      }
    }

    fn arb_length() -> impl Strategy<Value = f64> {
      proptest::num::f64::NORMAL | proptest::num::f64::SUBNORMAL | proptest::num::f64::ZERO
    }
  }
}
