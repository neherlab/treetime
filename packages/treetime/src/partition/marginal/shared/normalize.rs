use itertools::izip;
use ndarray::{Array, Array1, Array2, Axis, Dimension, Zip};
use treetime_primitives::LogLh;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm_inplace;
use treetime_utils::cast::ToF64Exact;

pub(crate) fn normalize_inplace(dis: &mut Array2<f64>) -> f64 {
  normalize_sites_inplace(dis).iter().sum()
}

pub(crate) fn normalize_sites_inplace(dis: &mut Array2<f64>) -> Array1<f64> {
  let norm = dis.sum_axis(Axis(1));
  let n_cols = dis.ncols().to_f64_exact();
  izip!(dis.outer_iter_mut(), &norm)
    .map(|(mut row, &n)| {
      if n > 0.0 && n.is_finite() {
        row /= n;
        n.ln()
      } else {
        row.fill(1.0 / n_cols);
        f64::NEG_INFINITY
      }
    })
    .collect()
}

pub(crate) fn normalize_sites_from_log(log_dis: Array2<f64>) -> (Array2<f64>, Array1<f64>) {
  let mut dis = log_dis;
  let site_log_norms = dis.rows_mut().into_iter().map(softmax_with_log_norm_inplace).collect();
  (dis, site_log_norms)
}

pub(crate) fn forward_log_lh_remove_child(node_log_lh: LogLh, child_log_lh: LogLh) -> LogLh {
  if node_log_lh == LogLh::IMPOSSIBLE && child_log_lh == LogLh::IMPOSSIBLE {
    LogLh::ZERO
  } else {
    LogLh::new(node_log_lh - child_log_lh)
  }
}

pub(crate) fn forward_log_lh_add_normalization(log_lh: LogLh, normalization: f64) -> LogLh {
  if normalization == f64::NEG_INFINITY {
    log_lh
  } else {
    log_lh + LogLh::new(normalization)
  }
}

pub(crate) fn divide_out<D: Dimension>(profile: &Array<f64, D>, message: &Array<f64, D>) -> Array<f64, D> {
  Zip::from(profile)
    .and(message)
    .map_collect(|&value, &divisor| divide_out_value(value, divisor))
}

pub(crate) fn divide_out_value(value: f64, divisor: f64) -> f64 {
  value / divisor.max(f64::MIN_POSITIVE)
}
