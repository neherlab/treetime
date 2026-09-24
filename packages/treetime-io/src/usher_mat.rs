use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use std::io::Write;
use std::path::Path;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub use util_usher_mat::{UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode};

pub fn usher_mat_pb_write_file(filepath: impl AsRef<Path>, tree: &UsherTree) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  usher_mat_pb_write(&mut f, tree)
    .wrap_err_with(|| format!("When writing Usher MAT protobuf file '{}'", filepath.display()))
}

fn usher_mat_pb_write(writer: &mut impl Write, tree: &UsherTree) -> Result<(), Report> {
  util_usher_mat::usher_mat_pb_write(writer, tree).wrap_err("When writing Usher MAT protobuf")
}

pub fn usher_mat_json_write_file(
  filepath: impl AsRef<Path>,
  tree: &UsherTree,
  options: &UsherMatJsonOptions,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  json_write_file(filepath, tree, JsonPretty(options.pretty))
    .wrap_err_with(|| format!("When writing Usher MAT JSON file: '{}'", filepath.display()))?;
  Ok(())
}

#[derive(SmartDefault)]
pub struct UsherMatJsonOptions {
  #[default = true]
  pretty: bool,
}
