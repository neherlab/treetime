use crate::error::OperationError;
use crate::make_internal_report;
use crate::make_report;
use statrs::distribution::{Discrete, Poisson};
use std::collections::BTreeMap;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::cast::usize_to_u64;

pub(crate) fn site_histogram(
  genome_length: usize,
  hits_per_site: &BTreeMap<usize, usize>,
) -> Result<SiteHistogram, OperationError> {
  if genome_length == 0 {
    return Err(OperationError::InvalidInput(make_report!(
      "The homoplasy statistics need at least one site, but the alignment is empty and no constant sites were given"
    )));
  }
  let unhit = genome_length.checked_sub(hits_per_site.len()).ok_or_else(|| {
    OperationError::InferenceFailed(make_internal_report!(
      "{} sites carry mutations, more than the {genome_length} sites of the genome",
      hits_per_site.len()
    ))
  })?;
  let max_hits = hits_per_site.values().copied().max().unwrap_or(0);
  let mut counts = vec![0_usize; max_hits + 1];
  counts[0] = unhit;
  for &hits in hits_per_site.values() {
    counts[hits] += 1;
  }
  let mutations: usize = hits_per_site.values().sum();
  let length = genome_length.to_f64_exact();

  if mutations == 0 {
    return Ok(SiteHistogram {
      genome_length,
      rows: vec![SiteHits {
        hits: 0,
        sites: genome_length,
        expected: length,
      }],
      log_likelihood_difference: 0.0,
    });
  }

  let poisson = Poisson::new(mutations.to_f64_exact() / length)
    .map_err(|err| OperationError::InferenceFailed(make_internal_report!("Poisson rate is invalid: {err}")))?;
  let ln_p = |k: usize| poisson.ln_pmf(usize_to_u64(k));

  let rows = counts
    .iter()
    .enumerate()
    .map(|(hits, &sites)| SiteHits {
      hits,
      sites,
      expected: length * ln_p(hits).exp(),
    })
    .collect();

  let observed: f64 = counts
    .iter()
    .enumerate()
    .map(|(hits, &sites)| sites.to_f64_exact() * ln_p(hits))
    .sum();
  let expected: f64 = (0..3 * counts.len())
    .map(|k| {
      let p = ln_p(k).exp();
      if p == 0.0 { 0.0 } else { p * ln_p(k) }
    })
    .sum();

  Ok(SiteHistogram {
    genome_length,
    rows,
    log_likelihood_difference: observed - length * expected,
  })
}

#[derive(Clone, Debug, PartialEq)]
pub struct SiteHistogram {
  pub genome_length: usize,
  pub rows: Vec<SiteHits>,
  pub log_likelihood_difference: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SiteHits {
  pub hits: usize,
  pub sites: usize,
  pub expected: f64,
}
