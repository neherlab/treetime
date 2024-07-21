use bytes::{Buf, BufMut};
use eyre::{Context, Report};
use prost::Message;
use std::io::{Read, Write};

pub mod mutation_detailed {
  include!(concat!(env!("OUT_DIR"), "/mutation_detailed.rs"));
}
pub mod parsimony {
  include!(concat!(env!("OUT_DIR"), "/parsimony.rs"));
}
pub mod taxodium {
  include!(concat!(env!("OUT_DIR"), "/taxodium.rs"));
}

pub fn usher_mat_pb_read_bytes(buf: impl Buf) -> Result<parsimony::Data, Report> {
  parsimony::Data::decode(buf).wrap_err("When decoding Usher MAT protobuf message")
}

pub fn usher_mat_pb_read(mut reader: impl Read) -> Result<parsimony::Data, Report> {
  let mut buf = Vec::new();
  reader.read_to_end(&mut buf)?;
  usher_mat_pb_read_bytes(&buf[..])
}

pub fn usher_mat_pb_write_bytes(buf: &mut impl BufMut, data: &parsimony::Data) -> Result<(), Report> {
  data.encode(buf).wrap_err("When encoding Usher MAT protobuf message")
}

pub fn usher_mat_pb_write(writer: &mut impl Write, data: &parsimony::Data) -> Result<(), Report> {
  let mut buf = Vec::<u8>::with_capacity(data.encoded_len());
  usher_mat_pb_write_bytes(&mut buf, data)?;
  writer.write_all(&buf).wrap_err("When writing encoded protobuf message")
}

pub fn usher_mat_json_read_str(s: impl AsRef<str>) -> Result<parsimony::Data, Report> {
  serde_json::from_str(s.as_ref()).wrap_err("When reading Usher MAT JSON string")
}

pub fn usher_mat_json_read(reader: impl Read) -> Result<parsimony::Data, Report> {
  serde_json::from_reader(reader).wrap_err("When reading Usher MAT JSON")
}

pub fn usher_mat_json_write_bytes(data: &parsimony::Data) -> Result<String, Report> {
  serde_json::to_string(data).wrap_err("When writing Usher MAT JSON string")
}

pub fn usher_mat_json_write(writer: &mut impl Write, data: &parsimony::Data) -> Result<(), Report> {
  serde_json::to_writer(writer, data).wrap_err("When writing Usher MAT JSON")
}
