mod convolution;
mod multiplication;
pub mod runner;

pub use convolution::ConvolutionRunner;
pub use multiplication::ChainMultiplicationRunner;
pub use multiplication::MultiplicationRunner;
pub use runner::TestRunner;
pub use runner::run_tests_generic;
