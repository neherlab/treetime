pub mod algorithm_summary;
pub mod algorithms;
pub mod console;
pub mod framework;
pub mod functions;
pub mod metrics;
pub mod output;

// Re-export commonly used items
pub use algorithm_summary::AlgorithmSummary;
pub use algorithms::ConvolutionAlgorithm;
pub use console::ConvolutionTestConsole;
pub use functions::exponential::{ExponentialFlatResult, ExponentialTestRunner};
pub use framework::{GenericConvolutionTestFramework, TestSummary};
pub use functions::gaussian::{GaussianFlatResult, GaussianTestRunner};
pub use output::{BaseFlatResult, TestOutputWriter, ToFlatResult};
