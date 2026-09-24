use bytes::BytesMut;
use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use std::io::Write;
use std::path::Path;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::json::{JsonPretty, json_write, json_write_file, json_write_str};

pub use util_usher_mat::{UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode};

pub fn usher_mat_pb_write_file(filepath: impl AsRef<Path>, tree: &UsherTree) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  usher_mat_pb_write(&mut f, tree)
    .wrap_err_with(|| format!("When writing Usher MAT protobuf file '{}'", filepath.display()))
}

pub(crate) fn usher_mat_pb_write_bytes(tree: &UsherTree) -> Result<Vec<u8>, Report> {
  let mut buf = BytesMut::new();
  util_usher_mat::usher_mat_pb_write_bytes(&mut buf, tree).wrap_err("When writing Usher MAT protobuf bytes")?;
  Ok(buf.to_vec())
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

pub fn usher_mat_json_write_str(tree: &UsherTree, options: &UsherMatJsonOptions) -> Result<String, Report> {
  json_write_str(tree, JsonPretty(options.pretty)).wrap_err("When writing Usher MAT JSON string")
}

pub fn usher_mat_json_write(
  writer: &mut impl Write,
  tree: &UsherTree,
  options: &UsherMatJsonOptions,
) -> Result<(), Report> {
  json_write(writer, tree, JsonPretty(options.pretty)).wrap_err("When writing Usher MAT JSON")
}

#[derive(SmartDefault)]
pub struct UsherMatJsonOptions {
  #[default = true]
  pretty: bool,
}
