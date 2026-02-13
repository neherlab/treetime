use ctor::ctor;
use eyre::Report;
use treetime_utils::init::global::global_init;
use treetime_validation::testing::run::run_validation_tests;

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  run_validation_tests()
}
