pub mod orchestrator;
pub mod results;
pub mod runner;
pub mod summary;
pub mod test_case;

pub use orchestrator::GenericConvolutionTestFramework;
pub use results::{TestFailure, TestResult, TestRunOutcome};
pub use runner::{ConvolutionTestRunner, TraitBasedTestRunner};
pub use summary::TestSummary;
pub use test_case::TestCase;
