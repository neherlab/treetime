pub mod types;

use crate::types::Phyloxml;
use quick_xml::DeError;
use quick_xml::de::from_reader;
use quick_xml::se::Serializer;
use serde::Serialize;
use std::io;

pub fn phyloxml_read(reader: impl io::Read) -> Result<Phyloxml, DeError> {
  let reader = io::BufReader::new(reader);
  from_reader(reader)
}

pub fn phyloxml_write(mut writer: impl io::Write, phyloxml: &Phyloxml) -> io::Result<()> {
  writer.write_all(b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")?;

  let mut document = String::new();
  let mut serializer = Serializer::with_root(&mut document, Some("phyloxml")).map_err(io::Error::other)?;
  serializer.indent(' ', 2);
  serializer.expand_empty_elements(true);
  phyloxml.serialize(serializer).map_err(io::Error::other)?;

  writer.write_all(document.as_bytes())
}

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;

  #[ctor(unsafe)]
  fn init() {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
