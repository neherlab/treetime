#[cfg(test)]
mod tests {
  use crate::usher_mat::usher_mat_pb_write_file;
  use bytes::BytesMut;
  use eyre::{Report, WrapErr};
  use pretty_assertions::assert_eq;
  use std::fs;
  use tempfile::TempDir;
  use util_usher_mat::{usher_mat_pb_read_bytes, usher_mat_pb_write_bytes};

  #[test]
  fn test_usher_mat_pb_write_file_matches_encoded_bytes() -> Result<(), Report> {
    let tree = helpers::tree();
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let path = dir.path().join("tree.pb");

    usher_mat_pb_write_file(&path, &tree)?;

    let written = fs::read(&path).wrap_err("When reading the written file")?;
    let mut expected = BytesMut::new();
    usher_mat_pb_write_bytes(&mut expected, &tree)?;
    assert_eq!(expected.to_vec(), written);
    Ok(())
  }

  #[test]
  fn test_usher_mat_pb_write_file_roundtrip() -> Result<(), Report> {
    let tree = helpers::tree();
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let path = dir.path().join("tree.pb");

    usher_mat_pb_write_file(&path, &tree)?;

    let written = fs::read(&path).wrap_err("When reading the written file")?;
    assert_eq!(tree, usher_mat_pb_read_bytes(&written[..])?);
    Ok(())
  }

  mod helpers {
    use crate::usher_mat::UsherTree;

    pub(super) fn tree() -> UsherTree {
      UsherTree {
        newick: "(A:1,B:2)root;".to_owned(),
        ..UsherTree::default()
      }
    }
  }
}
