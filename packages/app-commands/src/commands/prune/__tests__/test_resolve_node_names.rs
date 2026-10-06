#[cfg(test)]
mod tests {
  use crate::commands::prune::run::resolve_node_names;
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::o;

  #[test]
  fn test_resolve_node_names_prunes_the_first_node_of_a_repeated_name() {
    let tree = nwk_read(b"((A:0.1,B:0.1)X:0.1,(A:0.1,C:0.1)X:0.1)root;".as_slice()).unwrap();
    let names = tree.names();
    let keys_named = |name: &str| -> Vec<GraphNodeKey> {
      names
        .iter()
        .filter(|(_, node_name)| node_name.as_deref() == Some(name))
        .map(|(key, _)| *key)
        .sorted()
        .collect()
    };

    let actual = resolve_node_names(&btreeset! { o!("A"), o!("X"), o!("missing") }, &tree.graph, &names);

    assert_eq!((2, 2), (keys_named("A").len(), keys_named("X").len()));
    assert_eq!(btreeset! { keys_named("A")[0], keys_named("X")[0] }, actual);
  }
}
