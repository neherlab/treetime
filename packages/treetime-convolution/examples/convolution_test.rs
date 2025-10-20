use eyre::Report;
use treetime_convolution::testing::run::run_convolution_tests;

fn main() -> Result<(), Report> {
  run_convolution_tests()
}
