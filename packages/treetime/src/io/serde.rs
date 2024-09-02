use serde::Serializer;

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
