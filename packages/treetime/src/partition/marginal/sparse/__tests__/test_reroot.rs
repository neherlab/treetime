#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::marginal::sparse::reroot::reroot_sparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use crate::seq::mutation::Sub;
  use crate::test_utils::{deletion, find_edge_key, find_node_key_by_name, sparse_edge_obs};
  use eyre::Report;
  use helpers::c;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::reroot::{RerootResult, remove_stem_root};
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::Seq;

  #[test]
  fn test_sparse_reroot_stem_removal_gives_the_root_the_stem_child_sequence() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)R:0.001)STEM;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("named node exists");
    let (stem_key, r_key, a_key, b_key) = (key("STEM"), key("R"), key("A"), key("B"));
    let edge = |source: &str, target: &str| find_edge_key(&graph, &names, source, target).expect("edge exists");
    let (stem_edge, r_a_edge, r_b_edge) = (edge("STEM", "R"), edge("R", "A"), edge("R", "B"));
    let stem_seq = Seq::try_from_slice(b"ACGTACGT")?;
    let r_seq = Seq::try_from_slice(b"ACTTA--T")?;
    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet: alphabet.clone(),
      length: 8,
      root_sequence: stem_seq.clone(),
      obs_nodes: btreemap! {
        stem_key => SparseNodeObs::new(&stem_seq, &alphabet),
        r_key => SparseNodeObs::new(&r_seq, &alphabet),
        a_key => SparseNodeObs::new(&r_seq, &alphabet),
        b_key => SparseNodeObs::new(&r_seq, &alphabet),
      },
      obs_edges: btreemap! {
        stem_edge => sparse_edge_obs(
          vec![Sub::new(c(b'G'), 2_usize, c(b'T'))?],
          vec![deletion((5, 7), Seq::try_from_slice(b"CG")?)],
        ),
        r_a_edge => SparseEdgeObs::default(),
        r_b_edge => SparseEdgeObs::default(),
      },
    };
    let node_states = btreemap! {
      stem_key => SparseNodeState::leaf(&stem_seq),
      r_key => SparseNodeState::leaf(&r_seq),
      a_key => SparseNodeState::leaf(&r_seq),
      b_key => SparseNodeState::leaf(&r_seq),
    };
    let stem = remove_stem_root(&mut graph, stem_key)?.expect("STEM has one child");
    let changes = RerootResult {
      stem_removal: Some(stem),
      ..RerootResult::unchanged(r_key)
    };

    let recon = reroot_sparse(partition, jc69(JC69Params::default())?, node_states, &changes)?;

    assert_eq!(r_seq, recon.partition.root_sequence);
    assert_eq!(
      vec![r_key, a_key, b_key].into_iter().sorted().collect_vec(),
      recon.partition.obs_nodes.keys().copied().collect_vec()
    );
    assert_eq!(
      vec![r_a_edge, r_b_edge].into_iter().sorted().collect_vec(),
      recon.partition.obs_edges.keys().copied().collect_vec()
    );
    assert_eq!(
      vec![r_key, a_key, b_key].into_iter().sorted().collect_vec(),
      recon.node_states.keys().copied().collect_vec()
    );
    Ok(())
  }

  mod helpers {
    use treetime_primitives::AsciiChar;

    pub(super) fn c(byte: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(byte)
    }
  }
}
