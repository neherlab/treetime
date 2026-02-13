use indexmap::IndexMap;
use ndarray::Array1;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// Output empty value if false
///
/// Usage:
///     #[serde(serialize_with = "skip_serializing_if_false")]
//      pub is_foo: bool
#[allow(clippy::trivially_copy_pass_by_ref)]
pub fn skip_serializing_if_false<S>(value: &bool, serializer: S) -> Result<S::Ok, S::Error>
where
  S: Serializer,
{
  if *value {
    serializer.serialize_bool(*value)
  } else {
    serializer.serialize_none()
  }
}

/// Serialize Array1<T> as a simple JSON array
///
/// Usage:
///     #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
///     pub values: Array1<T>
pub fn array1_as_vec<T, S>(array: &Array1<T>, serializer: S) -> Result<S::Ok, S::Error>
where
  T: Serialize,
  S: Serializer,
{
  array.as_slice().unwrap().serialize(serializer)
}

/// Deserialize Array1<T> from a simple JSON array
///
/// Usage:
///     #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
///     pub values: Array1<T>
pub fn array1_from_vec<'de, T, D>(deserializer: D) -> Result<Array1<T>, D::Error>
where
  T: Deserialize<'de>,
  D: Deserializer<'de>,
{
  let vec = Vec::<T>::deserialize(deserializer)?;
  Ok(Array1::from_vec(vec))
}

/// Deserialize IndexMap<String, Array1<f64>> from JSON object
///
/// Usage:
///     #[serde(deserialize_with = "indexmap_array1_from_map")]
///     pub values: IndexMap<String, Array1<f64>>
pub fn indexmap_array1_from_map<'de, D>(deserializer: D) -> Result<IndexMap<String, Array1<f64>>, D::Error>
where
  D: Deserializer<'de>,
{
  let map = IndexMap::<String, Vec<f64>>::deserialize(deserializer)?;
  Ok(map.into_iter().map(|(k, v)| (k, Array1::from_vec(v))).collect())
}
