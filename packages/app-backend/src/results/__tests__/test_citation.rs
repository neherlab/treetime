#[cfg(test)]
mod tests {
  use crate::results::citation::citation;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_citation_links_its_doi() {
    let citation = citation();
    assert_eq!(
      ("10.1093/ve/vex042", "https://doi.org/10.1093/ve/vex042"),
      (citation.doi.as_str(), citation.url.as_str())
    );
  }
}
