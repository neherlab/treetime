use itertools::izip;
use ndarray::prelude::*;
use treetime_primitives::LogLh;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm;

pub fn normalize_inplace(dis: &mut Array2<f64>) -> f64 {
  let norm = dis.sum_axis(Axis(1));
  let n_cols = dis.ncols() as f64;
  let mut log_lh = 0.0;
  for (ri, mut row) in dis.outer_iter_mut().enumerate() {
    let n = norm[ri];
    if n > 0.0 && n.is_finite() {
      row /= n;
      log_lh += n.ln();
    } else {
      row.fill(1.0 / n_cols);
      log_lh += f64::NEG_INFINITY;
    }
  }
  log_lh
}

pub fn normalize_from_log(log_dis: &Array2<f64>) -> (Array2<f64>, f64) {
  let mut dis = Array2::zeros(log_dis.raw_dim());
  let mut total_log_lh = 0.0;

  for (mut out_row, log_row) in izip!(dis.rows_mut(), log_dis.rows()) {
    let (normalized, log_norm) = softmax_with_log_norm(log_row);
    out_row.assign(&normalized);
    total_log_lh += log_norm;
  }

  (dis, total_log_lh)
}

/// Remove a child's normalization scale from a forward cavity message.
///
/// A degenerate profile is represented by a uniform fallback with a
/// `f64::NEG_INFINITY` scale. Matching sentinels refer to the same removed
/// factor, so they cancel instead of performing the undefined IEEE-754
/// operation `-inf - (-inf)`.
pub(crate) fn forward_log_lh_remove_child(node_log_lh: LogLh, child_log_lh: LogLh) -> LogLh {
  if node_log_lh == LogLh::IMPOSSIBLE && child_log_lh == LogLh::IMPOSSIBLE {
    LogLh::ZERO
  } else {
    LogLh::new(node_log_lh - child_log_lh)
  }
}

/// Add a normalization scale to a forward cavity message.
///
/// The `f64::NEG_INFINITY` value marks a distribution replaced by the uniform
/// fallback. It has no finite scale to add to the conditional message. Other
/// non-finite values still propagate so unexpected NaN and positive infinity
/// remain observable.
pub(crate) fn forward_log_lh_add_normalization(log_lh: LogLh, normalization: f64) -> LogLh {
  if normalization == f64::NEG_INFINITY {
    log_lh
  } else {
    log_lh + LogLh::new(normalization)
  }
}
