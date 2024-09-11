use itertools::Itertools;
use serde::{Deserialize, Deserializer, Serializer};

/// Serde serializer for Letter sequences
pub fn serde_serialize_seq<S: Serializer>(seq: &[char], s: S) -> Result<S::Ok, S::Error> {
  let seq: String = seq.iter().collect();
  s.serialize_str(&seq)
}

/// Serde deserializer for Letter sequences
pub fn serde_deserialize_seq<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Vec<char>, D::Error> {
  let seq_str = String::deserialize(deserializer)?;
  let seq = seq_str.chars().collect_vec();
  Ok(seq)
}
