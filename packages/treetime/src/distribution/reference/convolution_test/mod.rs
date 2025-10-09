pub mod algorithms;
pub mod exponential;
pub mod framework;
pub mod gaussian;
pub mod output;

// Re-export commonly used items
pub use algorithms::ConvolutionAlgorithm;
pub use framework::{GenericConvolutionTestFramework, TestSummary};
pub use gaussian::{GaussianFlatResult, GaussianTestRunner};
pub use exponential::{ExponentialFlatResult, ExponentialTestRunner};
pub use output::{BaseFlatResult, TestOutputWriter, ToFlatResult};
