use crate::testing::metrics::config::DistributionConfig;
use crate::testing::metrics::distribution::histograms::{HistogramMetrics, compute_histogram_metrics};
use crate::testing::metrics::distribution::properties::{DistributionProperties, compute_distribution_properties};
use crate::testing::metrics::distribution::statistics::{StatisticalMetrics, compute_statistical_metrics};
use crate::testing::metrics::pointwise::errors::PointwiseErrors;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistributionMetrics {
  pub histograms: HistogramMetrics,
  pub statistics: StatisticalMetrics,
  pub properties: DistributionProperties,
}

impl DistributionMetrics {
  pub fn new(pointwise_errors: &PointwiseErrors) -> eyre::Result<Self> {
    Self::new_with_config(pointwise_errors, &DistributionConfig::default())
  }

  pub fn new_with_config(pointwise_errors: &PointwiseErrors, config: &DistributionConfig) -> eyre::Result<Self> {
    let histograms = compute_histogram_metrics(pointwise_errors, config)?;
    let statistics = compute_statistical_metrics(pointwise_errors)?;
    let properties = compute_distribution_properties(pointwise_errors)?;

    Ok(Self {
      histograms,
      statistics,
      properties,
    })
  }
}
