use crate::alphabet::alphabet::Alphabet;
use crate::build::insert::{Insertee, insert};
use crate::build::mat::Mat;
use crate::build::pipeline::encode_queries;
use crate::build::placement::place_exhaustive;
use crate::build::query::Query;
use crate::build::refresh::fitch_refresh;
use eyre::Report;
use itertools::Itertools;
use rand::seq::SliceRandom;
use std::collections::BTreeMap;
use treetime_primitives::Seq;
use treetime_utils::sync::random::get_random_number_generator;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RepruneParams {
  pub batches: usize,
  pub rounds: usize,
  pub seed: u64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RepruneStats {
  pub batches: usize,
  pub accepted: usize,
  pub gain: usize,
}

pub(crate) fn reprune(
  mat: Mat,
  seqs: &[Seq],
  alphabet: &Alphabet,
  params: &RepruneParams,
) -> Result<(Mat, RepruneStats), Report> {
  let mut rng = get_random_number_generator(params.seed);
  let mut stats = RepruneStats::default();
  let frame = mat.reference().clone();
  let encoded = encode_queries(seqs, &frame, alphabet);
  let mut mat = fitch_refresh(mat, &encoded, &frame, alphabet)?;
  for _ in 0..params.rounds {
    let mut samples = mat.samples().values().copied().collect_vec();
    samples.shuffle(&mut rng);
    let size = samples.len().div_ceil(params.batches.max(1)).max(1);
    for batch in samples.chunks(size) {
      stats.batches += 1;
      let before = mat.parsimony();
      let backup = mat.clone();
      let leaves: BTreeMap<usize, _> = mat.samples().iter().map(|(&leaf, &sample)| (sample, leaf)).collect();
      for sample in batch {
        mat.remove_sample(leaves[sample])?;
      }
      mat = fitch_refresh(mat, &encoded, &frame, alphabet)?;
      let queries = batch
        .iter()
        .map(|&sample| (sample, Query::from_seq(&seqs[sample], mat.reference(), alphabet)))
        .sorted_by_key(|(sample, query)| (query.n_missing(), *sample))
        .collect_vec();
      for (sample, query) in &queries {
        let placement = place_exhaustive(&mat, query, None, None, true, None)?;
        insert(&mut mat, query, placement.target, Insertee::Sample(*sample))?;
      }
      mat = fitch_refresh(mat, &encoded, &frame, alphabet)?;
      if mat.parsimony() < before {
        stats.accepted += 1;
        stats.gain += before - mat.parsimony();
      } else {
        mat = backup;
      }
    }
  }
  Ok((mat, stats))
}
