use indexmap::IndexMap;
use ndarray::{Array1, Array2};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

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

pub fn array1_as_vec<T, S>(array: &Array1<T>, serializer: S) -> Result<S::Ok, S::Error>
where
  T: Serialize,
  S: Serializer,
{
  serializer.collect_seq(array.iter())
}

pub fn array1_from_vec<'de, T, D>(deserializer: D) -> Result<Array1<T>, D::Error>
where
  T: Deserialize<'de>,
  D: Deserializer<'de>,
{
  let vec = Vec::<T>::deserialize(deserializer)?;
  Ok(Array1::from_vec(vec))
}

pub fn array2_as_vec<T, S>(array: &Array2<T>, serializer: S) -> Result<S::Ok, S::Error>
where
  T: Serialize,
  S: Serializer,
{
  let rows = array
    .rows()
    .into_iter()
    .map(|row| row.into_iter().collect::<Vec<_>>())
    .collect::<Vec<_>>();
  rows.serialize(serializer)
}

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

pub fn option_array1_as_vec<T, S>(array: &Option<Array1<T>>, serializer: S) -> Result<S::Ok, S::Error>
where
  T: Serialize,
  S: Serializer,
{
  match array {
    Some(arr) => serializer.collect_seq(arr.iter()),
    None => serializer.serialize_none(),
  }
}

pub fn option_array1_from_vec<'de, T, D>(deserializer: D) -> Result<Option<Array1<T>>, D::Error>
where
  T: Deserialize<'de>,
  D: Deserializer<'de>,
{
  Option::<Vec<T>>::deserialize(deserializer).map(|opt| opt.map(Array1::from_vec))
}

pub fn indexmap_array1_from_map<'de, D>(deserializer: D) -> Result<IndexMap<String, Array1<f64>>, D::Error>
where
  D: Deserializer<'de>,
{
  let map = IndexMap::<String, Vec<f64>>::deserialize(deserializer)?;
  Ok(map.into_iter().map(|(k, v)| (k, Array1::from_vec(v))).collect())
}
