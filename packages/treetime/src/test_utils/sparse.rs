use crate::partition::storage::sparse::SparseEdgeObs;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;

pub(crate) fn fitch_edge_obs(subs: Vec<Sub>) -> SparseEdgeObs {
  let mut obs = SparseEdgeObs::default();
  obs.set_fitch_subs(subs);
  obs
}

pub(crate) fn sparse_edge_obs(subs: Vec<Sub>, indels: Vec<InDel>) -> SparseEdgeObs {
  let mut obs = fitch_edge_obs(subs);
  obs.indels = indels;
  obs
}
