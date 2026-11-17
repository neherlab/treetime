use crate::gtr::gtr::{Mu, avg_transition};
use crate::progress::LogSink;
use crate::progress_warn;
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, Array2, Array3, ArrayView1, Axis};
use smart_default::SmartDefault;
use treetime_utils::array::ndarray::{is_max_above, outer};
use treetime_utils::cast::ToF64Exact;

pub const DEFAULT_PSEUDO_COUNT: f64 = 1.0;

const TINY_NUMBER: f64 = 1e-12;

pub(crate) fn infer_gtr_from_counts(
  counts: &MutationCounts,
  options: &InferGtrOptions,
  log: &dyn LogSink,
) -> Result<InferGtrResult, Report> {
  infer_gtr_group(&[counts], GtrSharing::default(), options, log)?
    .pop()
    .ok_or_else(|| make_internal_report!("The model fit of one partition returned no model"))
}

pub(crate) fn infer_gtr_group(
  counts: &[&MutationCounts],
  sharing: GtrSharing,
  options: &InferGtrOptions,
  log: &dyn LogSink,
) -> Result<Vec<InferGtrResult>, Report> {
  let InferGtrOptions {
    fixed_pi,
    pc,
    dp,
    max_iter,
  } = options;
  let Some(first) = counts.first() else {
    return make_internal_error!("A model fit needs the counts of at least one partition");
  };
  let N = first.Ti.len();
  if let Some(other) = counts.iter().find(|counts| counts.Ti.len() != N) {
    return make_internal_error!(
      "Partitions that share model parameters need one alphabet, but their state counts are {N} and {}",
      other.Ti.len()
    );
  }

  let pc_mat = {
    let mut pc_mat = Array2::from_elem((N, N), *pc);
    pc_mat.diag_mut().fill(0.0);
    pc_mat
  };
  let data = counts
    .iter()
    .map(|MutationCounts { nij, Ti, root_state }| {
      let mut nij = nij.clone();
      nij.diag_mut().fill(0.0);
      PartitionCounts { nij, Ti, root_state }
    })
    .collect_vec();
  let time_weights = time_weights(&data);

  let mut pi_old = vec![Array1::<f64>::zeros(N); data.len()];
  let mut pi = {
    let pi = fixed_pi.clone().unwrap_or_else(|| Array1::ones(N));
    let pi_sum = pi.sum();
    vec![pi / pi_sum; data.len()]
  };
  let mut W = vec![Array2::ones((N, N)); data.len()];
  let mut mu = data
    .iter()
    .map(|counts| (counts.nij.sum() + pc) / (counts.Ti.sum() + pc))
    .collect_vec();

  for _ in 0..*max_iter {
    if max_distance(&pi_old, &pi) < *dp {
      break;
    }
    pi_old.clone_from(&pi);

    W = if sharing.exchange_rates {
      let W = shared_exchange_rates(&data, &pi, &mu, &time_weights, &pc_mat)?;
      vec![W; data.len()]
    } else {
      izip!(&data, &pi, &mu)
        .map(|(counts, pi, &mu)| {
          let mut W = (&(&counts.nij.view() + &counts.nij.t() + 2.0 * &pc_mat) / mu)
            / (&outer(pi, counts.Ti)? + &outer(counts.Ti, pi)? + TINY_NUMBER + 2.0 * &pc_mat);
          W.diag_mut().fill(0.0);
          W /= avg_transition(&W, pi);
          Ok(W)
        })
        .collect::<Result<Vec<_>, Report>>()?
    };

    if fixed_pi.is_none() {
      let pseudo_counts = pc_mat.sum_axis(Axis(1));
      pi = if sharing.frequencies {
        let changes = data.iter().fold(Array1::<f64>::zeros(N), |sum, counts| {
          sum + counts.nij.sum_axis(Axis(1))
        });
        let roots = data
          .iter()
          .fold(Array1::<f64>::zeros(N), |sum, counts| sum + counts.root_state);
        let times = izip!(&data, &W, &mu).fold(Array1::<f64>::zeros(N), |sum, (counts, W, &mu)| {
          sum + mu * W.dot(counts.Ti)
        });
        let mut shared: Array1<f64> =
          (&changes + &pseudo_counts + &roots) / (TINY_NUMBER + times + roots.sum() + &pseudo_counts);
        shared /= shared.sum();
        vec![shared; data.len()]
      } else {
        izip!(&data, &W, &mu)
          .map(|(counts, W, &mu)| {
            let mut pi = (&counts.nij.sum_axis(Axis(1)) + &pseudo_counts + counts.root_state)
              / (TINY_NUMBER + mu * W.dot(counts.Ti) + counts.root_state.sum() + &pseudo_counts);
            pi /= pi.sum();
            pi
          })
          .collect()
      };
      mu = izip!(&data, &W, &pi)
        .map(|(counts, W, pi)| (counts.nij.sum() + pc) / (pi.dot(&W.dot(counts.Ti)) + pc))
        .collect();
    } else {
      mu = izip!(&data, &W, &pi)
        .map(|(counts, W, pi)| (counts.nij.sum() + pc) / (pi.dot(&W.dot(pi)) * counts.Ti.sum() + pc))
        .collect();
    }
  }

  if max_distance(&pi_old, &pi) > *dp {
    progress_warn!(
      log,
      "When inferring GTR parameters: The iterative scheme has not converged."
    );
  } else if pi.iter().any(|pi| (pi.sum() - 1.0).abs() > *dp) {
    progress_warn!(
      log,
      "When inferring GTR parameters: Proper normalization was not reached."
    );
  }
  Ok(
    izip!(W, pi, mu)
      .map(|(W, pi, mu)| InferGtrResult {
        W,
        pi,
        mu: Mu::FromMatrix(mu),
      })
      .collect(),
  )
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct GtrSharing {
  pub(crate) exchange_rates: bool,
  pub(crate) frequencies: bool,
}

struct PartitionCounts<'a> {
  nij: Array2<f64>,
  Ti: &'a Array1<f64>,
  root_state: &'a Array1<f64>,
}

fn time_weights(data: &[PartitionCounts<'_>]) -> Vec<f64> {
  let times = data.iter().map(|counts| counts.Ti.sum()).collect_vec();
  let total = times.iter().sum::<f64>();
  if total > 0.0 {
    times.iter().map(|time| time / total).collect()
  } else {
    vec![1.0 / data.len().to_f64_exact(); data.len()]
  }
}

fn shared_exchange_rates(
  data: &[PartitionCounts<'_>],
  pi: &[Array1<f64>],
  mu: &[f64],
  time_weights: &[f64],
  pc_mat: &Array2<f64>,
) -> Result<Array2<f64>, Report> {
  let shape = pc_mat.raw_dim();
  let changes = data.iter().fold(Array2::<f64>::zeros(shape), |sum, counts| {
    sum + &counts.nij + counts.nij.t()
  });
  let times = izip!(data, pi, mu).try_fold(Array2::<f64>::zeros(shape), |sum, (counts, pi, &mu)| {
    Ok::<_, Report>(sum + mu * (&outer(pi, counts.Ti)? + &outer(counts.Ti, pi)?))
  })?;
  let mean_mu = izip!(time_weights, mu).map(|(weight, mu)| weight * mu).sum::<f64>();
  let mut W = (&changes + 2.0 * pc_mat) / (&times + mean_mu * (TINY_NUMBER + 2.0 * pc_mat));
  W.diag_mut().fill(0.0);
  let mean_rate = izip!(time_weights, pi)
    .map(|(weight, pi)| weight * avg_transition(&W, pi))
    .sum::<f64>();
  W /= mean_rate;
  Ok(W)
}

fn max_distance(pi_old: &[Array1<f64>], pi: &[Array1<f64>]) -> f64 {
  izip!(pi_old, pi)
    .map(|(pi_old, pi)| distance(pi_old, pi))
    .fold(0.0, f64::max)
}

#[derive(Clone, Debug)]
pub struct MutationCounts {
  pub(crate) nij: Array2<f64>,

  pub(crate) Ti: Array1<f64>,

  pub(crate) root_state: Array1<f64>,
}

#[derive(Clone, Debug, PartialEq, SmartDefault)]
pub struct InferGtrOptions {
  pub(crate) fixed_pi: Option<Array1<f64>>,

  #[default(DEFAULT_PSEUDO_COUNT)]
  pub(crate) pc: f64,

  #[default = 1e-5]
  pub(crate) dp: f64,

  #[default = 40]
  pub(crate) max_iter: usize,
}

#[derive(Clone, Debug)]
pub struct InferGtrResult {
  pub(crate) W: Array2<f64>,
  pub(crate) pi: Array1<f64>,
  pub(crate) mu: Mu,
}

pub(crate) fn distance(pi_old: &Array1<f64>, pi: &Array1<f64>) -> f64 {
  (pi_old - pi).mapv(|x| x * x).sum().sqrt()
}

#[expect(
  clippy::needless_pass_by_value,
  reason = "array views are Copy; clippy misses the Copy impl of ndarray 0.17 views"
)]
pub(crate) fn is_profile_informative(profile: ArrayView1<'_, f64>, n_states: usize) -> bool {
  let uniform_threshold = 1.0 / n_states.to_f64_exact() + 1e-10;
  is_max_above(&profile, uniform_threshold)
}

pub(crate) fn get_branch_mutation_matrix(
  msg_to_child: &Array2<f64>,
  msg_to_parent: &Array2<f64>,
  exp_qt: &Array2<f64>,
) -> Array3<f64> {
  let (n_sites, n_states) = msg_to_parent.dim();
  let mut result = Array3::zeros((n_sites, n_states, n_states));

  for a in 0..n_sites {
    let mut site_sum = 0.0;

    for i in 0..n_states {
      for j in 0..n_states {
        let val = msg_to_parent[[a, i]] * exp_qt[[i, j]] * msg_to_child[[a, j]];
        result[[a, i, j]] = val;
        site_sum += val;
      }
    }

    if site_sum > 0.0 {
      for i in 0..n_states {
        for j in 0..n_states {
          result[[a, i, j]] /= site_sum;
        }
      }
    }
  }

  result
}

pub(crate) fn accumulate_mutation_counts(
  mut_stack: &Array3<f64>,
  branch_length: f64,
  nij: &mut Array2<f64>,
  Ti: &mut Array1<f64>,
) {
  let (n_sites, n_states, _) = mut_stack.dim();

  for i in 0..n_states {
    for j in 0..n_states {
      let mut sum = 0.0;
      for a in 0..n_sites {
        sum += mut_stack[[a, i, j]];
      }
      nij[[i, j]] += sum;
    }
  }

  for k in 0..n_states {
    let mut parent_sum = 0.0;
    let mut child_sum = 0.0;
    for a in 0..n_sites {
      for i in 0..n_states {
        parent_sum += mut_stack[[a, i, k]];
      }
      for j in 0..n_states {
        child_sum += mut_stack[[a, k, j]];
      }
    }
    Ti[k] += 0.5 * branch_length * (parent_sum + child_sum);
  }
}
