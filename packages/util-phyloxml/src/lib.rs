pub mod types;

use crate::types::Phyloxml;
use deser::Error;
use deser_xml::{Indent, SerializerConfig};
use std::io;

const DOCUMENT: SerializerConfig = SerializerConfig::new().declaration(true).indent(Indent::Spaces(2));

pub fn phyloxml_read(reader: impl io::Read) -> Result<Phyloxml, Error> {
  deser_xml::from_reader(reader)
}

pub fn phyloxml_write(writer: impl io::Write, phyloxml: &Phyloxml) -> Result<(), Error> {
  DOCUMENT.to_writer(writer, phyloxml)
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
