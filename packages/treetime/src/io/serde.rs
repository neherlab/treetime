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

/// Serialize Array1<f64> as a simple JSON array
///
/// Usage:
///     #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
///     pub values: Array1<f64>
pub fn array1_as_vec<S>(array: &Array1<f64>, serializer: S) -> Result<S::Ok, S::Error>
where
  S: Serializer,
{
  array.as_slice().unwrap().serialize(serializer)
}

/// Deserialize Array1<f64> from a simple JSON array
///
/// Usage:
///     #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
///     pub values: Array1<f64>
pub fn array1_from_vec<'de, D>(deserializer: D) -> Result<Array1<f64>, D::Error>
where
  D: Deserializer<'de>,
{
  let vec = Vec::<f64>::deserialize(deserializer)?;
  Ok(Array1::from_vec(vec))
}
