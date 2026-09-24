use crate::partition::storage::sparse::SparseEdgeObs;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;

pub(crate) fn sparse_edge_obs(subs: Vec<Sub>, indels: Vec<InDel>) -> SparseEdgeObs {
  let mut obs = SparseEdgeObs::with_fitch_subs(subs);
  obs.indels = indels;
  obs
}
