use crate::clock::clock_regression::ClockParams;
use crate::make_report;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime_io::fasta::FastaRecord;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum TimeMarginalMode {
  #[default]
  Never,
  Always,
  OnlyFinal,
}

/// Compute effective time_marginal mode, promoting Never to OnlyFinal when
/// confidence is requested with rate uncertainty prerequisites.
///
/// v0 (wrappers.py:478):
///   time_marginal = 'confidence-only' if (calc_confidence and time_marginal == 'never') else time_marginal
pub fn compute_effective_time_marginal(
  time_marginal: TimeMarginalMode,
  confidence: bool,
  clock_std_dev: Option<f64>,
  covariation: bool,
) -> TimeMarginalMode {
  if confidence && time_marginal == TimeMarginalMode::Never {
    if clock_std_dev.is_some() || covariation {
      info!("--confidence: promoting time-marginal from never to only-final for CI estimation");
      TimeMarginalMode::OnlyFinal
    } else {
      warn!(
        "Cannot estimate confidence intervals without clock rate uncertainty. \
         Specify --clock-std-dev or rerun with --covariation. \
         Proceeding without confidence estimation."
      );
      TimeMarginalMode::Never
    }
  } else {
    time_marginal
  }
}

/// Build covariation-aware ClockParams when covariation is enabled.
///
/// v0 (clock_tree.py:277-285):
///   branch_variance = (max(0, clock_length) + tip_slack^2 * om) * om   [leaves]
///   branch_variance = max(0, clock_length) * om                         [internal]
/// where om = 1/seq_len, tip_slack = OVER_DISPERSION = 10 (config.py:8)
///
/// Mapping to v1 ClockParams:
///   variance_factor = 1/seq_len, variance_offset = 0, variance_offset_leaf = tip_slack^2/seq_len^2
pub fn build_covariation_clock_params(
  covariation: bool,
  sequence_length: Option<usize>,
  tip_slack: Option<f64>,
  aln: Option<&[FastaRecord]>,
) -> Result<Option<ClockParams>, Report> {
  if !covariation {
    return Ok(None);
  }

  let seq_len = if let Some(aln_data) = aln {
    get_common_length(aln_data)? as f64
  } else {
    sequence_length.ok_or_else(|| make_report!("--sequence-length required for --covariation without alignment"))?
      as f64
  };

  let tip_slack = tip_slack.unwrap_or(10.0);

  info!("Covariation-aware clock regression: seq_len={seq_len}, tip_slack={tip_slack}");

  Ok(Some(ClockParams {
    variance_factor: 1.0 / seq_len,
    variance_offset: 0.0,
    variance_offset_leaf: tip_slack * tip_slack / (seq_len * seq_len),
  }))
}
