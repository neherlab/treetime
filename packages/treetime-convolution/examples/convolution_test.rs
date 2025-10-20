use ctor::ctor;
use eyre::Report;
use treetime_convolution::testing::run::run_convolution_tests;
use treetime_utils::global_init::global_init;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  run_convolution_tests()
}
