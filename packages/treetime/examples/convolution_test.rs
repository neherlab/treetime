use eyre::Report;
use treetime::distribution::reference::convolution_test::run::run_convolution_tests;

fn main() -> Result<(), Report> {
  run_convolution_tests()
}
