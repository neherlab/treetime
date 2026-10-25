#[cfg(test)]
mod tests {
  use self::generators::{Op, gen_ops};
  use crate::key_map::NodeMap;
  use crate::node::GraphNodeKey;
  use itertools::Itertools;
  use proptest::{prop_assert_eq, proptest};
  use std::collections::BTreeMap;

  proptest! {
    #[test]
    fn test_prop_key_map_matches_btree_map(ops in gen_ops()) {
      let mut map = NodeMap::new();
      let mut oracle = BTreeMap::new();
      let mut answers = (vec![], vec![]);
      for op in ops {
        match op {
          Op::Insert(index, value) => {
            answers.0.push(map.insert(GraphNodeKey(index), value));
            answers.1.push(oracle.insert(GraphNodeKey(index), value));
          },
          Op::Remove(index) => {
            answers.0.push(map.remove(GraphNodeKey(index)));
            answers.1.push(oracle.remove(&GraphNodeKey(index)));
          },
          Op::RetainBelow(limit) => {
            map.retain(|_, value| *value < limit);
            oracle.retain(|_, value| *value < limit);
          },
        }
      }

      let entries = map.iter().map(|(key, &value)| (key, value)).collect_vec();
      let expected = oracle.iter().map(|(&key, &value)| (key, value)).collect_vec();
      prop_assert_eq!(&answers.1, &answers.0);
      prop_assert_eq!((expected.len(), expected), (map.len(), entries));
    }
  }

  mod generators {
    use proptest::collection::vec;
    use proptest::prelude::{Just, Strategy, prop_oneof};

    #[derive(Clone, Debug)]
    pub(super) enum Op {
      Insert(usize, u8),
      Remove(usize),
      RetainBelow(u8),
    }

    pub(super) fn gen_ops() -> impl Strategy<Value = Vec<Op>> {
      let op = prop_oneof![
        4 => (0_usize..40, 0_u8..100).prop_map(|(index, value)| Op::Insert(index, value)),
        2 => (0_usize..40).prop_map(Op::Remove),
        1 => (0_u8..100).prop_map(Op::RetainBelow),
        1 => Just(Op::RetainBelow(u8::MAX)),
      ];
      vec(op, 0..60)
    }
  }
}
