#[cfg(test)]
mod tests {
  use crate::{Phyloxml, phyloxml_read, phyloxml_write};
  use indoc::indoc;
  use pretty_assertions::assert_eq;

  fn write_to_string(doc: &Phyloxml) -> String {
    let mut buf = Vec::new();
    phyloxml_write(&mut buf, doc).unwrap();
    String::from_utf8(buf).unwrap()
  }

  #[test]
  fn test_phyloxml_write_emits_declaration_indent_and_order() {
    let input = r#"<phyloxml><phylogeny rooted="true"><name>test tree</name></phylogeny></phyloxml>"#;
    let doc = phyloxml_read(input.as_bytes()).unwrap();
    let actual = write_to_string(&doc);
    let expected = indoc! {r#"
      <?xml version="1.0" encoding="UTF-8"?>
      <phyloxml>
        <phylogeny rooted="true">
          <name>test tree</name>
        </phylogeny>
      </phyloxml>"#};
    assert_eq!(expected, actual);
  }
}
