use eyre::{Report, WrapErr};
use std::io::{self, Write};
use treetime_utils::init::openblas::get_openblas_info_str;

pub(crate) fn print_debug_info() -> Result<(), Report> {
  writeln!(io::stdout().lock(), "{}", get_openblas_info_str())
    .wrap_err("When writing debug information to standard output")
}
