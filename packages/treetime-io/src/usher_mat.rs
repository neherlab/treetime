use eyre::{Report, WrapErr};
use std::path::Path;
use treetime_utils::io::file::write_file_with;
use util_usher_mat::usher_mat_pb_write;
pub use util_usher_mat::{UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode};

pub fn usher_mat_pb_write_file(filepath: impl AsRef<Path>, tree: &UsherTree) -> Result<(), Report> {
  write_file_with(filepath, |writer| {
    usher_mat_pb_write(writer, tree).wrap_err("When writing UShER MAT protobuf")
  })
}
