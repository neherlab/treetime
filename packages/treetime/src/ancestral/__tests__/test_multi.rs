use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::multi::{MarginalPartitionParams, PartitionPlan, reconstruct_marginal_partitions};
use crate::ancestral::sample::SampleMode;
use crate::gtr::get_gtr::GtrModelName;
use crate::payload::ancestral::GraphAncestral;
use pretty_assertions::assert_eq;
use treetime_io::nwk::nwk_read_str;

/// Two amino-acid partitions of different lengths reconstruct on one shared tree in a single
/// multi-partition traversal. Each is independent (its own length), and the stop codon `*` is carried
/// as a real state rather than rejected (the bug when reconstruction used the 20-state no-stop
/// alphabet). Fully in-memory: the graph is parsed from a string and sequences are built directly.
#[test]
fn test_multi_reconstructs_each_cds_independently_with_stop_codon() {
  let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.1)root;").unwrap();
  let aa = Alphabet::new(AlphabetName::Aa).unwrap();

  let plans = vec![
    helpers::plan("S", &aa, &[("A", "MC*"), ("B", "MA*")]),
    helpers::plan("N", &aa, &[("A", "KL"), ("B", "KM")]),
  ];
  let params = MarginalPartitionParams {
    dense: Some(false),
    reconstruct_tip_states: false,
    sample_from_profile: SampleMode::default(),
    seed: None,
    ignore_missing_alns: false,
  };

  let reconstructed = reconstruct_marginal_partitions(&graph, plans, &params).unwrap();

  let names = reconstructed.iter().map(|p| p.name.clone()).collect::<Vec<_>>();
  assert_eq!(vec!["S".to_owned(), "N".to_owned()], names);

  let s_root = reconstructed[0].partition.read_arc().root_sequence(&graph).unwrap();
  assert_eq!(3, s_root.len());
  assert_eq!(Some('*'), s_root.iter().last().map(|c| char::from(*c)));

  let n_root = reconstructed[1].partition.read_arc().root_sequence(&graph).unwrap();
  assert_eq!(2, n_root.len());
  assert_eq!(Some('K'), n_root.iter().next().map(|c| char::from(*c)));
}

mod helpers {
  use super::*;
  use treetime_io::fasta::FastaRecord;
  use treetime_primitives::Seq;

  pub fn plan(name: &str, alphabet: &Alphabet, seqs: &[(&str, &str)]) -> PartitionPlan {
    PartitionPlan {
      name: name.to_owned(),
      alphabet: alphabet.clone(),
      gtr_model: GtrModelName::Infer,
      sequences: seqs
        .iter()
        .enumerate()
        .map(|(index, (seq_name, seq))| FastaRecord {
          seq_name: (*seq_name).to_owned(),
          desc: None,
          seq: Seq::try_from_str(seq).unwrap(),
          index,
        })
        .collect(),
      annotation: None,
      reference_override: None,
    }
  }
}
