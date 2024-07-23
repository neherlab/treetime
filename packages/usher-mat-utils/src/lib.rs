use bytes::{Buf, BufMut};
use eyre::{Context, Report};
use prost::Message;
use std::io::{Read, Write};

mod mutation_detailed {
  include!(concat!(env!("OUT_DIR"), "/mutation_detailed.rs"));
}
mod parsimony {
  include!(concat!(env!("OUT_DIR"), "/parsimony.rs"));
}
mod taxodium {
  include!(concat!(env!("OUT_DIR"), "/taxodium.rs"));
}

pub type UsherTree = parsimony::Data;
pub type UsherTreeNode = parsimony::CondensedNode;
pub type UsherMutation = parsimony::Mut;
pub type UsherMutationList = parsimony::MutationList;
pub type UsherMetadata = parsimony::NodeMetadata;

pub fn usher_mat_pb_read_bytes(buf: impl Buf) -> Result<UsherTree, Report> {
  UsherTree::decode(buf).wrap_err("When decoding Usher MAT protobuf message")
}

pub fn usher_mat_pb_read(mut reader: impl Read) -> Result<UsherTree, Report> {
  let mut buf = Vec::new();
  reader.read_to_end(&mut buf)?;
  usher_mat_pb_read_bytes(&buf[..])
}

pub fn usher_mat_pb_write_bytes(buf: &mut impl BufMut, tree: &UsherTree) -> Result<(), Report> {
  tree.encode(buf).wrap_err("When encoding Usher MAT protobuf message")
}

pub fn usher_mat_pb_write(writer: &mut impl Write, tree: &UsherTree) -> Result<(), Report> {
  let mut buf = Vec::<u8>::with_capacity(tree.encoded_len());
  usher_mat_pb_write_bytes(&mut buf, tree)?;
  writer.write_all(&buf).wrap_err("When writing encoded protobuf message")
}

pub fn usher_mat_json_read_str(s: impl AsRef<str>) -> Result<UsherTree, Report> {
  serde_json::from_str(s.as_ref()).wrap_err("When reading Usher MAT JSON string")
}

pub fn usher_mat_json_read(reader: impl Read) -> Result<UsherTree, Report> {
  serde_json::from_reader(reader).wrap_err("When reading Usher MAT JSON")
}

pub fn usher_mat_json_write_str(tree: &UsherTree) -> Result<String, Report> {
  serde_json::to_string(tree).wrap_err("When writing Usher MAT JSON string")
}

pub fn usher_mat_json_write(writer: &mut impl Write, tree: &UsherTree) -> Result<(), Report> {
  serde_json::to_writer(writer, tree).wrap_err("When writing Usher MAT JSON")
}
