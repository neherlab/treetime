#[cfg(test)]
mod tests {
  use crate::__tests__::test_graph_support::tests::helpers::names_on_star;
  use crate::node::GraphNodeKey;
  use crate::pair_by_name::{NamePairing, key_by_name, pair_by_name};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_pair_by_name_keeps_the_first_entry_of_a_repeated_name() -> Result<(), Report> {
    let (a, b) = (GraphNodeKey(0), GraphNodeKey(1));
    let names = names_on_star(btreemap! { a => Some(o!("A")), b => Some(o!("B")) })?;
    let entries = vec![(o!("A"), 1), (o!("B"), 2), (o!("A"), 3), (o!("B"), 4), (o!("A"), 5)];

    let actual = pair_by_name([a, b], &names, entries);

    assert_eq!(
      NamePairing {
        by_node: btreemap! { a => 1, b => 2 },
        unmatched: vec![],
        duplicate_entry_names: vec_of_owned!["A", "B"],
      },
      actual
    );
    Ok(())
  }

  #[test]
  fn test_pair_by_name_lists_unmatched_entries_in_input_order() -> Result<(), Report> {
    let a = GraphNodeKey(0);
    let names = names_on_star(btreemap! { a => Some(o!("A")) })?;
    let entries = vec![(o!("Z"), 1), (o!("A"), 2), (o!("C"), 3), (o!("Z"), 4)];

    let actual = pair_by_name([a], &names, entries);

    assert_eq!(
      NamePairing {
        by_node: btreemap! { a => 2 },
        unmatched: vec![(o!("Z"), 1), (o!("C"), 3)],
        duplicate_entry_names: vec_of_owned!["Z"],
      },
      actual
    );
    Ok(())
  }

  #[test]
  fn test_pair_by_name_leaves_out_candidates_without_an_entry() -> Result<(), Report> {
    let (a, b, unnamed) = (GraphNodeKey(0), GraphNodeKey(1), GraphNodeKey(2));
    let names = names_on_star(btreemap! { a => Some(o!("A")), b => Some(o!("B")), unnamed => None })?;

    let actual = pair_by_name([a, b, unnamed], &names, vec![(o!("A"), 1)]);

    assert_eq!(
      NamePairing {
        by_node: btreemap! { a => 1 },
        unmatched: vec![],
        duplicate_entry_names: vec![],
      },
      actual
    );
    Ok(())
  }

  #[test]
  fn test_pair_by_name_gives_the_entry_to_every_candidate_with_its_name() -> Result<(), Report> {
    let (a1, b, a2) = (GraphNodeKey(0), GraphNodeKey(1), GraphNodeKey(2));
    let names = names_on_star(btreemap! { a1 => Some(o!("A")), b => Some(o!("B")), a2 => Some(o!("A")) })?;

    let actual = pair_by_name([a1, b, a2], &names, vec![(o!("A"), 1), (o!("B"), 2)]);

    assert_eq!(
      NamePairing {
        by_node: btreemap! { a1 => 1, b => 2, a2 => 1 },
        unmatched: vec![],
        duplicate_entry_names: vec![],
      },
      actual
    );
    Ok(())
  }

  #[test]
  fn test_pair_by_name_ignores_nodes_that_are_not_candidates() -> Result<(), Report> {
    let (leaf, internal) = (GraphNodeKey(0), GraphNodeKey(1));
    let names = names_on_star(btreemap! { leaf => Some(o!("A")), internal => Some(o!("I")) })?;

    let actual = pair_by_name([leaf], &names, vec![(o!("I"), 1)]);

    assert_eq!(
      NamePairing {
        by_node: btreemap! {},
        unmatched: vec![(o!("I"), 1)],
        duplicate_entry_names: vec![],
      },
      actual
    );
    Ok(())
  }

  #[test]
  fn test_pair_by_name_key_by_name_resolves_the_written_name_of_each_node() -> Result<(), Report> {
    let (b, a1, a2) = (GraphNodeKey(0), GraphNodeKey(1), GraphNodeKey(2));
    let names = names_on_star(btreemap! { b => Some(o!("B")), a1 => Some(o!("A")), a2 => Some(o!("A")) })?;

    let actual = key_by_name([a2, b, a1], &names);

    assert_eq!(btreemap! { "A" => a1, "A_2" => a2, "B" => b }, actual);
    Ok(())
  }
}
