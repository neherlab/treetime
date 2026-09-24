#![allow(
  clippy::derive_partial_eq_without_eq,
  reason = "Prost messages derive PartialEq without Eq"
)]

use eyre::{Context, Report};
use prost::Message;
use prost::bytes::{Buf, BufMut};
use std::collections::BTreeSet;
use std::io::{Read, Write};

mod mutation_detailed {
  #![allow(
    dead_code,
    unnameable_types,
    unreachable_pub,
    unknown_lints,
    no_comments,
    topological_ordering,
    clippy::all,
    reason = "Prost-generated module"
  )]
  include!(concat!(env!("OUT_DIR"), "/mutation_detailed.rs"));
}
mod parsimony {
  #![allow(
    dead_code,
    unnameable_types,
    unreachable_pub,
    unknown_lints,
    no_comments,
    topological_ordering,
    clippy::all,
    reason = "Prost-generated module"
  )]
  include!(concat!(env!("OUT_DIR"), "/parsimony.rs"));
}
mod taxodium {
  #![allow(
    dead_code,
    unnameable_types,
    unreachable_pub,
    unknown_lints,
    no_comments,
    topological_ordering,
    clippy::all,
    reason = "Prost-generated module"
  )]
  include!(concat!(env!("OUT_DIR"), "/taxodium.rs"));
}

pub type UsherTreeNode = parsimony::CondensedNode;
pub type UsherMutation = parsimony::Mut;
pub type UsherMutationList = parsimony::MutationList;
pub type UsherMetadata = parsimony::NodeMetadata;
pub fn usher_mat_pb_read(mut reader: impl Read) -> Result<UsherTree, Report> {
  let mut buf = Vec::new();
  reader
    .read_to_end(&mut buf)
    .wrap_err("When reading Usher MAT protobuf input")?;
  usher_mat_pb_read_bytes(&buf[..])
}

pub fn usher_mat_pb_read_bytes(buf: impl Buf) -> Result<UsherTree, Report> {
  UsherTree::decode(buf).wrap_err("When decoding Usher MAT protobuf message")
}

pub fn usher_mat_pb_write(writer: &mut impl Write, tree: &UsherTree) -> Result<(), Report> {
  let mut buf = Vec::<u8>::with_capacity(tree.encoded_len());
  usher_mat_pb_write_bytes(&mut buf, tree)?;
  writer.write_all(&buf).wrap_err("When writing encoded protobuf message")
}

pub fn usher_mat_pb_write_bytes(buf: &mut impl BufMut, tree: &UsherTree) -> Result<(), Report> {
  tree.encode(buf).wrap_err("When encoding Usher MAT protobuf message")
}

pub type UsherTree = parsimony::Data;

impl UsherTree {
  pub fn get_all_positions(&self) -> BTreeSet<i32> {
    self
      .node_mutations
      .iter()
      .flat_map(|ml| ml.mutation.iter().map(|m| m.position))
      .collect()
  }
}
