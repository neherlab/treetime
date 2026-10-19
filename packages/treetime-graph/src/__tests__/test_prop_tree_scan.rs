#[cfg(test)]
mod tests {
  use self::generators::gen_parents;
  use self::helpers::{record_scan_down, record_scan_up, sorted_keys, visits_before_children, visits_before_parent};
  use super::super::test_graph_support::tests::graph_from_parents;
  use crate::node::GraphNodeKey;
  use crate::tree_view::TreeView;
  use itertools::Itertools;
  use proptest::prelude::*;
  use std::collections::BTreeMap;

  proptest! {
    #[test]
    fn test_prop_tree_scan_up_visits_each_node_once_after_its_children_in_tree_order(
      parents in gen_parents(),
    ) {
      let graph = graph_from_parents(&parents).unwrap();
      let view = TreeView::new(&graph).unwrap();

      let single = record_scan_up(&view, 1).unwrap();
      let multi = record_scan_up(&view, 4).unwrap();

      let expected = view
        .preorder()
        .iter()
        .map(|&key| (key, view.children(key).iter().map(|&(child, _)| child).collect_vec()))
        .collect::<BTreeMap<_, _>>();
      prop_assert_eq!(&expected, &single.received);
      prop_assert_eq!(&single.received, &multi.received);
      prop_assert_eq!(sorted_keys(&view), multi.visits.iter().copied().sorted().collect_vec());
      prop_assert_eq!(Vec::<GraphNodeKey>::new(), visits_before_children(&view, &multi));
      prop_assert_eq!((1, 1), (single.live_messages, multi.live_messages));
    }

    #[test]
    fn test_prop_tree_scan_down_visits_each_node_once_after_its_parent(
      parents in gen_parents(),
    ) {
      let graph = graph_from_parents(&parents).unwrap();
      let view = TreeView::new(&graph).unwrap();

      let single = record_scan_down(&view, 1).unwrap();
      let multi = record_scan_down(&view, 4).unwrap();

      let expected = view
        .preorder()
        .iter()
        .map(|&key| (key, view.parent(key).map(|(parent, _)| parent).into_iter().collect_vec()))
        .collect::<BTreeMap<_, _>>();
      prop_assert_eq!(&expected, &single.received);
      prop_assert_eq!(&single.received, &multi.received);
      prop_assert_eq!(sorted_keys(&view), multi.visits.iter().copied().sorted().collect_vec());
      prop_assert_eq!(Vec::<GraphNodeKey>::new(), visits_before_parent(&view, &multi));
      prop_assert_eq!((1, 1), (single.live_messages, multi.live_messages));
    }
  }

  mod generators {
    use itertools::Itertools;
    use proptest::prelude::*;

    pub(super) fn gen_parents() -> impl Strategy<Value = Vec<usize>> {
      (1_usize..60).prop_flat_map(|count| (0..count).map(|index| 0..index.max(1)).collect_vec())
    }
  }

  mod helpers {
    use crate::node::GraphNodeKey;
    use crate::tree_view::TreeView;
    use eyre::Report;
    use itertools::Itertools;
    use parking_lot::Mutex;
    use rayon::ThreadPoolBuilder;
    use std::collections::BTreeMap;
    use std::sync::Arc;

    pub(super) struct ScanRecord {
      pub(super) received: BTreeMap<GraphNodeKey, Vec<GraphNodeKey>>,
      pub(super) visits: Vec<GraphNodeKey>,
      pub(super) live_messages: usize,
    }

    pub(super) fn record_scan_up(view: &TreeView<'_>, threads: usize) -> Result<ScanRecord, Report> {
      let token = Arc::new(());
      let visits = Mutex::new(vec![]);
      let inputs = view.preorder().iter().map(|&key| (key, ()));
      let received = ThreadPoolBuilder::new().num_threads(threads).build()?.install(|| {
        view.scan_up(inputs, |key, (), children: Vec<Message>| {
          visits.lock().push(key);
          let received = children.iter().map(|message| message.key).collect_vec();
          Ok((received, Message::new(key, &token)))
        })
      })?;
      Ok(ScanRecord {
        received: received.into_map(),
        visits: visits.into_inner(),
        live_messages: Arc::strong_count(&token),
      })
    }

    pub(super) fn record_scan_down(view: &TreeView<'_>, threads: usize) -> Result<ScanRecord, Report> {
      let token = Arc::new(());
      let visits = Mutex::new(vec![]);
      let inputs = view.preorder().iter().map(|&key| (key, ()));
      let received = ThreadPoolBuilder::new().num_threads(threads).build()?.install(|| {
        view.scan_down(inputs, |key, (), parent: Option<Message>| {
          visits.lock().push(key);
          let received = parent.iter().map(|message| message.key).collect_vec();
          Ok((received, Message::new(key, &token)))
        })
      })?;
      Ok(ScanRecord {
        received: received.into_map(),
        visits: visits.into_inner(),
        live_messages: Arc::strong_count(&token),
      })
    }

    pub(super) fn sorted_keys(view: &TreeView<'_>) -> Vec<GraphNodeKey> {
      view.preorder().iter().copied().sorted().collect_vec()
    }

    pub(super) fn visits_before_children(view: &TreeView<'_>, record: &ScanRecord) -> Vec<GraphNodeKey> {
      let order = visit_positions(record);
      view
        .preorder()
        .iter()
        .copied()
        .filter(|key| view.children(*key).iter().any(|(child, _)| order[child] > order[key]))
        .collect_vec()
    }

    pub(super) fn visits_before_parent(view: &TreeView<'_>, record: &ScanRecord) -> Vec<GraphNodeKey> {
      let order = visit_positions(record);
      view
        .preorder()
        .iter()
        .copied()
        .filter(|key| view.parent(*key).is_some_and(|(parent, _)| order[&parent] > order[key]))
        .collect_vec()
    }

    fn visit_positions(record: &ScanRecord) -> BTreeMap<GraphNodeKey, usize> {
      record
        .visits
        .iter()
        .enumerate()
        .map(|(position, &key)| (key, position))
        .collect()
    }

    #[derive(Clone)]
    struct Message {
      key: GraphNodeKey,
      _token: Arc<()>,
    }

    impl Message {
      fn new(key: GraphNodeKey, token: &Arc<()>) -> Self {
        Self {
          key,
          _token: Arc::clone(token),
        }
      }
    }
  }
}
