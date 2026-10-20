#[cfg(test)]
mod tests {
  use self::helpers::{inputs_after_receiver, record_scan_down, record_scan_up, sorted_keys};
  use super::super::test_graph_support::tests::generators::gen_dag;
  use super::super::test_graph_support::tests::helpers::graph_from_dag;
  use crate::node::GraphNodeKey;
  use crate::schedule::GraphSchedule;
  use itertools::Itertools;
  use proptest::{prop_assert_eq, proptest};
  use std::collections::BTreeMap;

  proptest! {
    #[test]
    fn test_prop_schedule_scan_up_visits_each_node_once_after_its_children_in_edge_order(dag in gen_dag(60)) {
      let graph = graph_from_dag(&dag).unwrap();
      let schedule = GraphSchedule::new(&graph).unwrap();

      let single = record_scan_up(&schedule, 1).unwrap();
      let multi = record_scan_up(&schedule, 4).unwrap();

      let expected = schedule
        .preorder()
        .iter()
        .map(|&key| (key, schedule.children(key).iter().map(|&(child, _)| child).collect_vec()))
        .collect::<BTreeMap<_, _>>();
      prop_assert_eq!(&expected, &single.received);
      prop_assert_eq!(&single.received, &multi.received);
      prop_assert_eq!(sorted_keys(&schedule), multi.visits.iter().copied().sorted().collect_vec());
      prop_assert_eq!(Vec::<GraphNodeKey>::new(), inputs_after_receiver(&multi));
      prop_assert_eq!((1, 1), (single.live_messages, multi.live_messages));
    }

    #[test]
    fn test_prop_schedule_scan_down_visits_each_node_once_after_its_parents_in_edge_order(dag in gen_dag(60)) {
      let graph = graph_from_dag(&dag).unwrap();
      let schedule = GraphSchedule::new(&graph).unwrap();

      let single = record_scan_down(&schedule, 1).unwrap();
      let multi = record_scan_down(&schedule, 4).unwrap();

      let expected = schedule
        .preorder()
        .iter()
        .map(|&key| (key, schedule.parents(key).iter().map(|&(parent, _)| parent).collect_vec()))
        .collect::<BTreeMap<_, _>>();
      prop_assert_eq!(&expected, &single.received);
      prop_assert_eq!(&single.received, &multi.received);
      prop_assert_eq!(sorted_keys(&schedule), multi.visits.iter().copied().sorted().collect_vec());
      prop_assert_eq!(Vec::<GraphNodeKey>::new(), inputs_after_receiver(&multi));
      prop_assert_eq!((1, 1), (single.live_messages, multi.live_messages));
    }
  }

  mod helpers {
    use crate::node::GraphNodeKey;
    use crate::schedule::GraphSchedule;
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

    pub(super) fn record_scan_up(schedule: &GraphSchedule<'_>, threads: usize) -> Result<ScanRecord, Report> {
      let token = Arc::new(());
      let visits = Mutex::new(vec![]);
      let inputs = schedule.preorder().iter().map(|&key| (key, ()));
      let received = ThreadPoolBuilder::new().num_threads(threads).build()?.install(|| {
        schedule.scan_up(inputs, |key, (), children: Vec<Message>| {
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

    pub(super) fn record_scan_down(schedule: &GraphSchedule<'_>, threads: usize) -> Result<ScanRecord, Report> {
      let token = Arc::new(());
      let visits = Mutex::new(vec![]);
      let inputs = schedule.preorder().iter().map(|&key| (key, ()));
      let received = ThreadPoolBuilder::new().num_threads(threads).build()?.install(|| {
        schedule.scan_down(inputs, |key, (), parents: Vec<Message>| {
          visits.lock().push(key);
          let received = parents.iter().map(|message| message.key).collect_vec();
          Ok((received, Message::new(key, &token)))
        })
      })?;
      Ok(ScanRecord {
        received: received.into_map(),
        visits: visits.into_inner(),
        live_messages: Arc::strong_count(&token),
      })
    }

    pub(super) fn sorted_keys(schedule: &GraphSchedule<'_>) -> Vec<GraphNodeKey> {
      schedule.preorder().iter().copied().sorted().collect_vec()
    }

    pub(super) fn inputs_after_receiver(record: &ScanRecord) -> Vec<GraphNodeKey> {
      let order = record
        .visits
        .iter()
        .enumerate()
        .map(|(position, &key)| (key, position))
        .collect::<BTreeMap<_, _>>();
      record
        .received
        .iter()
        .filter(|(key, senders)| senders.iter().any(|sender| order[sender] > order[*key]))
        .map(|(&key, _)| key)
        .collect_vec()
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
