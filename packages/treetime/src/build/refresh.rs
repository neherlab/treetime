use crate::alphabet::alphabet::Alphabet;
use crate::build::mat::Mat;
use crate::make_internal_report;
use crate::partition::fitch::passes::create_fitch_partition;
use crate::seq::alignment::NodeSeqInput;
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::izip;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::Seq;

pub(crate) fn fitch_refresh(mat: Mat, seqs: &[Seq], alphabet: &Alphabet) -> Result<Mat, Report> {
  let root = mat.root();
  let old_reference = mat.reference().clone();
  let (mut graph, _, samples) = mat.into_parts();
  graph.build()?;
  let inputs = samples
    .iter()
    .map(|(&leaf, &sample)| {
      let seq = seqs[sample]
        .as_slice()
        .iter()
        .map(|&state| {
          if alphabet.is_gap(state) {
            alphabet.unknown()
          } else {
            state
          }
        })
        .collect();
      (
        leaf,
        NodeSeqInput {
          name: sample.to_string(),
          seq: Some(seq),
        },
      )
    })
    .collect::<BTreeMap<_, _>>();
  let fitch = create_fitch_partition(&graph, alphabet.clone(), inputs)?;
  let schedule = GraphSchedule::new(&graph)?;
  let mut subs = vec![
    vec![];
    graph
      .get_edges()
      .map(|edge| edge.key().as_usize() + 1)
      .max()
      .unwrap_or(0)
  ];
  let mut pending: BTreeMap<GraphNodeKey, (Seq, usize)> = BTreeMap::new();
  let mut reference = old_reference.clone();
  for &node in schedule.preorder() {
    let fitch_seq = &fitch
      .nodes
      .get(&node)
      .ok_or_else(|| make_internal_report!("Fitch reconstruction has no node {node}"))?
      .seq
      .sequence;
    let parent = schedule.parents(node).first().copied();
    let parent_seq = match parent {
      Some((parent, _)) => &pending[&parent].0,
      None => &old_reference,
    };
    let seq: Seq = izip!(fitch_seq.as_slice(), parent_seq.as_slice())
      .map(|(&state, &inherited)| if alphabet.is_canonical(state) { state } else { inherited })
      .collect();
    if let Some((parent, edge)) = parent {
      subs[edge.as_usize()] = izip!(parent_seq.as_slice(), seq.as_slice())
        .enumerate()
        .filter(|(_, (from, to))| from != to)
        .map(|(pos, (&from, &to))| Sub::new(from, pos, to))
        .collect::<Result<Vec<_>, Report>>()?;
      let remaining = &mut pending
        .get_mut(&parent)
        .ok_or_else(|| make_internal_report!("Parent {parent} is not pending"))?
        .1;
      *remaining -= 1;
      if *remaining == 0 {
        pending.remove(&parent);
      }
    } else {
      reference = seq.clone();
    }
    let n_children = schedule.children(node).len();
    if n_children > 0 {
      pending.insert(node, (seq, n_children));
    }
  }
  Ok(Mat::from_parts(graph, root, reference, subs, samples))
}
