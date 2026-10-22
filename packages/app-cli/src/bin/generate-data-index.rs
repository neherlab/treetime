use app_cli::data_index::generate_data_index;
use ctor::ctor;
use eyre::{Report, WrapErr};
use std::env;
use std::io::{Read, stdin};
use std::path::PathBuf;
use treetime_utils::init::global::global_init;
use treetime_utils::io::json::{JsonPretty, json_write_file};

#[ctor(unsafe)]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  let mut args = env::args().skip(1);
  let data_dir = args.next().map_or_else(|| PathBuf::from("data"), PathBuf::from);
  let output = args.next().map_or_else(|| PathBuf::from("-"), PathBuf::from);
  let mut file_list = String::new();
  stdin()
    .read_to_string(&mut file_list)
    .wrap_err("When reading the file list from standard input")?;
  let index = generate_data_index(&data_dir, &file_list)?;
  json_write_file(&output, &index, JsonPretty(true))
}
