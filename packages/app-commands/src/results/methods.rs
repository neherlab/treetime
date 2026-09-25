use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
use crate::results::timetree::{CoalescentPrior, TimetreeEstimates};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use treetime_utils::datetime::year_fraction::year_fraction_to_datestring;

const CITATION_TEXT: &str = "Sagulenko P, Puller V, Neher RA. TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution 4 (2018), vex042.";

const CITATION_DOI: &str = "https://doi.org/10.1093/ve/vex042";

/// The publication to cite for TreeTime.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct Citation {
  /// Reference in text form.
  pub text: String,
  /// Link to the publication.
  pub doi: String,
}

pub fn citation() -> Citation {
  Citation {
    text: CITATION_TEXT.to_owned(),
    doi: CITATION_DOI.to_owned(),
  }
}

pub fn timetree_methods(version: &str, config: &TreetimeTimetreeArgsRaw, estimates: &TimetreeEstimates) -> String {
  [
    format!(
      "A time-scaled phylogeny of {} samples was inferred with TreeTime {version} (timetree command).",
      estimates.samples
    ),
    rate_sentence(config, estimates),
    filter_sentence(config, estimates),
    prior_sentence(estimates.coalescent_prior),
    estimates.relaxed_clock.map_or_else(String::new, |relax| {
      format!(
        "A relaxed clock was used (slack {}, coupling {}).",
        relax.slack, relax.coupling
      )
    }),
    root_sentence(estimates),
    format!("Please cite: {CITATION_TEXT}"),
  ]
  .into_iter()
  .filter(|sentence| !sentence.is_empty())
  .collect::<Vec<_>>()
  .join(" ")
}

fn rate_sentence(config: &TreetimeTimetreeArgsRaw, estimates: &TimetreeEstimates) -> String {
  let Some(rate) = estimates.clock_rate else {
    return String::new();
  };
  if estimates.clock_rate_fixed {
    let spread = config.clock_std_dev.map_or_else(String::new, |std| {
      format!(" with standard deviation {}", rate_text(std))
    });
    format!(
      "The clock rate was fixed at {} substitutions per site per year{spread}.",
      rate_text(rate)
    )
  } else {
    let spread = estimates
      .clock_rate_std
      .map_or_else(String::new, |std| format!(" (standard deviation {})", rate_text(std)));
    format!(
      "The clock rate was estimated at {} substitutions per site per year{spread}.",
      rate_text(rate)
    )
  }
}

fn filter_sentence(config: &TreetimeTimetreeArgsRaw, estimates: &TimetreeEstimates) -> String {
  let filter = if config.clock_filter > 0.0 {
    format!(
      "Samples whose root-to-tip residual exceeded {} interquartile distances were treated as clock outliers. ",
      config.clock_filter
    )
  } else {
    String::new()
  };
  format!(
    "{filter}{} of {} samples had no usable date or were clock outliers and did not constrain the clock model.",
    estimates.excluded_samples, estimates.samples
  )
}

fn prior_sentence(prior: CoalescentPrior) -> String {
  let described = match prior {
    CoalescentPrior::None => return String::new(),
    CoalescentPrior::Fixed { tc } => format!("constant size, Tc = {tc} years"),
    CoalescentPrior::Optimized => "constant size, optimized Tc".to_owned(),
    CoalescentPrior::Skyline { points, stiffness } => format!("skyline, {points} points, stiffness {stiffness}"),
  };
  format!("A coalescent prior was used ({described}).")
}

fn root_sentence(estimates: &TimetreeEstimates) -> String {
  let Some(date) = estimates.root_date else {
    return String::new();
  };
  let range = estimates.root_interval.map_or_else(String::new, |interval| {
    format!(
      " (90% interval {} to {})",
      year_fraction_to_datestring(interval.lower),
      year_fraction_to_datestring(interval.upper)
    )
  });
  format!("The root was dated to {}{range}.", year_fraction_to_datestring(date))
}

fn rate_text(rate: f64) -> String {
  format!("{rate:.2e}")
}
