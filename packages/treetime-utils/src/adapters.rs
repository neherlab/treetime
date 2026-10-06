use deser::adapters::{DeserializeAs, FromInto, SerializeAs, TryFromInto};
use deser::de::SinkHandle;
use deser::ser::{Chunk, SeqEmitter, SerializeHandle};
use deser::{Atom, ContainerShape, Deserialize, Error, Serialize, State};
use ndarray::iter::Iter;
use ndarray::{Array1, Array2, ArrayView1, Ix1, ShapeError};

pub struct ArrayVec;

impl<T: Serialize> SerializeAs<Array1<T>> for ArrayVec {
  fn serialize_as<'a>(value: &'a Array1<T>, state: &mut State) -> Result<Chunk<'a>, Error> {
    Ok(Chunk::seq(Elements(value.iter()), state))
  }

  fn container_shape_as(value: &Array1<T>) -> ContainerShape {
    ContainerShape::new().with_len(value.len())
  }
}

impl<'de, T: Deserialize<'de> + Send + 'static> DeserializeAs<'de, Array1<T>> for ArrayVec {
  fn deserialize_into_as<'out>(out: &'out mut Option<Array1<T>>, state: &mut State) -> SinkHandle<'out, 'de> {
    <FromInto<Vec<T>> as DeserializeAs<'de, Array1<T>>>::deserialize_into_as(out, state)
  }
}

pub struct Array2Rows;

impl<T: Serialize> SerializeAs<Array2<T>> for Array2Rows {
  fn serialize_as<'a>(value: &'a Array2<T>, state: &mut State) -> Result<Chunk<'a>, Error> {
    Ok(Chunk::seq(RowsEmitter { array: value, row: 0 }, state))
  }

  fn container_shape_as(value: &Array2<T>) -> ContainerShape {
    ContainerShape::new().with_len(value.nrows())
  }
}

impl<'de, T: Deserialize<'de> + Send + 'static> DeserializeAs<'de, Array2<T>> for Array2Rows {
  fn deserialize_into_as<'out>(out: &'out mut Option<Array2<T>>, state: &mut State) -> SinkHandle<'out, 'de> {
    <TryFromInto<Rows<T>> as DeserializeAs<'de, Array2<T>>>::deserialize_into_as(out, state)
  }
}

pub struct TrueOrNull;

impl SerializeAs<bool> for TrueOrNull {
  fn serialize_as<'a>(value: &'a bool, _state: &mut State) -> Result<Chunk<'a>, Error> {
    Ok(Chunk::Atom(if *value { Atom::Bool(true) } else { Atom::Null }))
  }
}

impl<'de> DeserializeAs<'de, bool> for TrueOrNull {
  fn deserialize_into_as<'out>(out: &'out mut Option<bool>, state: &mut State) -> SinkHandle<'out, 'de> {
    <FromInto<MaybeTrue> as DeserializeAs<'de, bool>>::deserialize_into_as(out, state)
  }

  fn initial_value_as() -> Option<bool> {
    Some(false)
  }
}

struct Elements<'a, T>(Iter<'a, T, Ix1>);

impl<T: Serialize> SeqEmitter for Elements<'_, T> {
  fn next(&mut self, _state: &mut State) -> Result<Option<SerializeHandle<'_>>, Error> {
    Ok(self.0.next().map(SerializeHandle::to))
  }
}

struct RowsEmitter<'a, T> {
  array: &'a Array2<T>,
  row: usize,
}

impl<T: Serialize> SeqEmitter for RowsEmitter<'_, T> {
  fn next(&mut self, state: &mut State) -> Result<Option<SerializeHandle<'_>>, Error> {
    if self.row == self.array.nrows() {
      return Ok(None);
    }
    let row = Row(self.array.row(self.row));
    self.row += 1;
    Ok(Some(SerializeHandle::arena(row, state)))
  }
}

struct Row<'a, T>(ArrayView1<'a, T>);

impl<T: Serialize> Serialize for Row<'_, T> {
  fn serialize(&self, state: &mut State) -> Result<Chunk<'_>, Error> {
    Ok(Chunk::seq(Elements(self.0.iter()), state))
  }

  fn container_shape(&self) -> ContainerShape {
    ContainerShape::new().with_len(self.0.len())
  }
}

#[derive(Deserialize)]
struct Rows<T>(Vec<Vec<T>>);

impl<T> TryFrom<Rows<T>> for Array2<T> {
  type Error = ShapeError;

  fn try_from(Rows(rows): Rows<T>) -> Result<Self, ShapeError> {
    let shape = (rows.len(), rows.first().map_or(0, Vec::len));
    Array2::from_shape_vec(shape, rows.into_iter().flatten().collect())
  }
}

#[derive(Deserialize)]
struct MaybeTrue(Option<bool>);

impl From<MaybeTrue> for bool {
  fn from(MaybeTrue(value): MaybeTrue) -> Self {
    value.unwrap_or(false)
  }
}
