pub mod algorithm_summary;
pub mod algorithms;
pub mod console;
pub mod exponential;
pub mod framework;
pub mod gaussian;
pub mod metrics;
pub mod output;

// Re-export commonly used items
pub use algorithm_summary::AlgorithmSummary;
pub use algorithms::ConvolutionAlgorithm;
pub use console::ConvolutionTestConsole;
pub use exponential::{ExponentialFlatResult, ExponentialTestRunner};
pub use framework::{GenericConvolutionTestFramework, TestSummary};
pub use gaussian::{GaussianFlatResult, GaussianTestRunner};
pub use output::{BaseFlatResult, TestOutputWriter, ToFlatResult};
