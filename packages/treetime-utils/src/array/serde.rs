use indexmap::IndexMap;
use ndarray::{Array1, Array2};
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

/// Serialize Array2<T> as a nested JSON array (row-major)
///
/// Usage:
///     #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
///     pub values: Array2<T>
pub fn array2_as_vec<T, S>(array: &Array2<T>, serializer: S) -> Result<S::Ok, S::Error>
where
  T: Serialize,
  S: Serializer,
{
  let rows: Vec<&[T]> = array.rows().into_iter().map(|row| row.to_slice().unwrap()).collect();
  rows.serialize(serializer)
}

/// Deserialize Array2<T> from a nested JSON array (row-major)
///
/// Usage:
///     #[serde(serialize_with = "array2_as_vec", deserialize_with = "array2_from_vec")]
///     pub values: Array2<T>
pub fn array2_from_vec<'de, T, D>(deserializer: D) -> Result<Array2<T>, D::Error>
where
  T: Deserialize<'de>,
  D: Deserializer<'de>,
{
  let nested = Vec::<Vec<T>>::deserialize(deserializer)?;
  let nrows = nested.len();
  let ncols = nested.first().map_or(0, Vec::len);
  let flat: Vec<T> = nested.into_iter().flatten().collect();
  Array2::from_shape_vec((nrows, ncols), flat).map_err(serde::de::Error::custom)
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
