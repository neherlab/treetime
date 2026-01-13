pub mod aggregate;
pub mod config;
pub mod distribution;
pub mod metrics;
pub mod pointwise;
pub mod spatial;

pub use metrics::ConvolutionMetrics;

pub type AccuracyMetrics = ConvolutionMetrics;
