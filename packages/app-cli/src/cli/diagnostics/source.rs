#![allow(
  clippy::wildcard_enum_match_arm,
  reason = "application layer: counts and indices to f64, integer division for averaging, graph node access and CLI setup invariants, and default variant matches"
)]

use color_eyre::Section;
use eyre::Report;
use itertools::Itertools;
use miette::{Diagnostic, LabeledSpan, NamedSource, Severity, SourceCode, SourceSpan};
use saphyr::{LoadableYamlNode, MarkedYaml, Scalar, YamlData};
use serde_json::Value;
use serde_saphyr::DuplicateKeyPolicy;
use std::collections::BTreeMap;
use std::fmt::{self, Display, Formatter};
use treetime_utils::{make_error, make_report};

/// A config file's text plus a JSON-pointer to source-span index for caret placement.
///
/// The span index is built from a `saphyr` parse of the same text, so it maps every node (and each
/// mapping key) to a byte range in the original document. Diagnostics carry a JSON pointer; the
/// renderer looks the pointer up here to draw a caret. When the text cannot be parsed for spans, the
/// index is empty and diagnostics render without carets rather than failing.
pub struct ConfigSource {
  name: String,
  text: String,
  spans: BTreeMap<String, NodeSpan>,
}

impl ConfigSource {
  /// Build a source and its span index from a filename and the raw document text.
  pub fn new(name: impl Into<String>, text: impl Into<String>) -> Self {
    let name = name.into();
    let text = text.into();
    let mut spans = BTreeMap::new();
    if let Ok(docs) = MarkedYaml::load_from_str(&text) {
      if let Some(doc) = docs.first() {
        let table = char_byte_table(&text);
        index_node(&mut spans, &table, doc, "");
      }
    }
    Self { name, text, spans }
  }

  /// Span of the value at `pointer`, if known.
  pub fn span_for(&self, pointer: &str) -> Option<SourceSpan> {
    self.spans.get(pointer).map(|node| node.value)
  }

  /// Span of the mapping key that introduces `pointer`, falling back to the value span.
  pub fn key_span_for(&self, pointer: &str) -> Option<SourceSpan> {
    self.spans.get(pointer).map(|node| node.key.unwrap_or(node.value))
  }

  /// A miette source over the whole document, for snippet rendering.
  pub fn named_source(&self) -> NamedSource<String> {
    NamedSource::new(&self.name, self.text.clone())
  }
}

/// A pending diagnostic, before its JSON pointer is resolved to a concrete span.
///
/// Passes collect these; `render_and_bail` turns each into a rendered `ConfigDiagnostic`. The pointer
/// is resolved late so the same pass logic works whether or not spans are available.
pub struct RawDiagnostic {
  pub pointer: Option<String>,
  pub use_key_span: bool,
  pub code: String,
  pub message: String,
  pub label: String,
  pub help: Option<String>,
}

impl RawDiagnostic {
  /// A diagnostic with a diagnostic code and headline message; attach location and help fluently.
  #[must_use]
  pub fn new(code: impl Into<String>, message: impl Into<String>) -> Self {
    Self {
      pointer: None,
      use_key_span: false,
      code: code.into(),
      message: message.into(),
      label: "here".to_owned(),
      help: None,
    }
  }

  /// Point the diagnostic at a JSON pointer into the document (its value span).
  #[must_use]
  pub fn at(mut self, pointer: impl Into<String>) -> Self {
    self.pointer = Some(pointer.into());
    self
  }

  /// Draw the caret under the mapping key rather than the value (for unknown/misused keys).
  #[must_use]
  pub fn key_span(mut self) -> Self {
    self.use_key_span = true;
    self
  }

  /// Attach an actionable `help:` line (suggestions, valid values).
  #[must_use]
  pub fn help(mut self, help: impl Into<String>) -> Self {
    self.help = Some(help.into());
    self
  }

  /// Resolve the pointer to a span against `source` and build the renderable diagnostic.
  fn resolve(self, source: &ConfigSource) -> ConfigDiagnostic {
    let span = self.pointer.as_deref().and_then(|pointer| {
      if self.use_key_span {
        source.key_span_for(pointer)
      } else {
        source.span_for(pointer)
      }
    });
    ConfigDiagnostic {
      message: self.message,
      code: self.code,
      help: self.help,
      src: source.named_source(),
      span,
      label: self.label,
    }
  }
}

/// Build an eyre error for a batch of config diagnostics.
///
/// The error chain is a terse, stable headline, `"{top_message}: {problems}"`, so callers and tests
/// can assert on it directly. The full caret-annotated source report rides along as a color-eyre
/// section, so the globally installed handler prints it once, with location and backtrace intact,
/// rather than the diagnostics writing to stderr on their own. An empty diagnostic list is success.
pub fn render_and_bail(source: &ConfigSource, top_message: &str, diags: Vec<RawDiagnostic>) -> Result<(), Report> {
  if diags.is_empty() {
    return Ok(());
  }
  let related: Vec<ConfigDiagnostic> = diags.into_iter().map(|diag| diag.resolve(source)).collect();
  let problems = related.iter().map(|diag| diag.message.as_str()).join("; ");
  let report = ConfigReport {
    message: top_message.to_owned(),
    related,
  };

  let mut rendered = String::new();
  miette::GraphicalReportHandler::new()
    .render_report(&mut rendered, &report)
    .map_err(|err| make_report!("could not render config diagnostics: {err}"))?;

  make_error!("{top_message}: {problems}").with_section(move || rendered.trim_end().to_owned())
}

/// Parse a config document (JSON or YAML), returning an error carrying a syntax diagnostic on failure.
///
/// YAML is a superset of JSON, so one parser reads both. Parsing rejects duplicate mapping keys and
/// non-finite floats (`.inf`, `.nan`): a configuration carries neither, so each is a hard parse error.
/// A parse failure becomes a caret-annotated report against `source`, attached to the returned error,
/// matching how every other config problem is surfaced.
pub fn parse_config_document(source: &ConfigSource, text: &str) -> Result<Value, Report> {
  let options = serde_saphyr::options! {
    duplicate_keys: DuplicateKeyPolicy::Error,
    reject_non_finite_typeless_float: true,
  };
  match serde_saphyr::from_str_with_options::<Value>(text, options) {
    Ok(value) => Ok(value),
    Err(err) => {
      render_and_bail(
        source,
        "invalid configuration",
        vec![RawDiagnostic::new(
          "config::syntax",
          format!("could not parse config: {err}"),
        )],
      )?;
      unreachable!("render_and_bail returns an error whenever diagnostics are present");
    },
  }
}

/// One rendered diagnostic: a message, code, optional caret, and optional help, over the document.
#[derive(Debug)]
struct ConfigDiagnostic {
  message: String,
  code: String,
  help: Option<String>,
  src: NamedSource<String>,
  span: Option<SourceSpan>,
  label: String,
}

impl Display for ConfigDiagnostic {
  fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
    f.write_str(&self.message)
  }
}

impl std::error::Error for ConfigDiagnostic {}

impl Diagnostic for ConfigDiagnostic {
  fn code<'a>(&'a self) -> Option<Box<dyn Display + 'a>> {
    Some(Box::new(self.code.clone()))
  }

  fn severity(&self) -> Option<Severity> {
    Some(Severity::Error)
  }

  fn help<'a>(&'a self) -> Option<Box<dyn Display + 'a>> {
    let help = self.help.clone()?;
    Some(Box::new(help))
  }

  fn source_code(&self) -> Option<&dyn SourceCode> {
    Some(&self.src)
  }

  fn labels(&self) -> Option<Box<dyn Iterator<Item = LabeledSpan> + '_>> {
    let span = self.span?;
    let label = LabeledSpan::new_with_span(Some(self.label.clone()), span);
    Some(Box::new(std::iter::once(label)))
  }
}

/// The batched report: a headline plus every collected diagnostic as a related entry.
#[derive(Debug)]
struct ConfigReport {
  message: String,
  related: Vec<ConfigDiagnostic>,
}

impl Display for ConfigReport {
  fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
    f.write_str(&self.message)
  }
}

impl std::error::Error for ConfigReport {}

impl Diagnostic for ConfigReport {
  fn severity(&self) -> Option<Severity> {
    Some(Severity::Error)
  }

  fn related<'a>(&'a self) -> Option<Box<dyn Iterator<Item = &'a dyn Diagnostic> + 'a>> {
    Some(Box::new(self.related.iter().map(|diag| -> &dyn Diagnostic { diag })))
  }
}

/// One node's spans: the value range, and the key range when the node is a mapping member.
struct NodeSpan {
  value: SourceSpan,
  key: Option<SourceSpan>,
}

/// Walk a `saphyr` node tree, recording a span for every node keyed by its JSON pointer.
fn index_node(spans: &mut BTreeMap<String, NodeSpan>, table: &[usize], node: &MarkedYaml, pointer: &str) {
  spans.insert(
    pointer.to_owned(),
    NodeSpan {
      value: span_of(table, node),
      key: None,
    },
  );

  match &node.data {
    YamlData::Mapping(map) => {
      for (key_node, value_node) in map {
        let Some(key) = node_key(key_node) else {
          continue;
        };
        let child = format!("{pointer}/{}", escape_pointer(&key));
        index_node(spans, table, value_node, &child);
        if let Some(node) = spans.get_mut(&child) {
          node.key = Some(span_of(table, key_node));
        }
      }
    },
    YamlData::Sequence(items) => {
      for (position, item) in items.iter().enumerate() {
        index_node(spans, table, item, &format!("{pointer}/{position}"));
      }
    },
    _ => {},
  }
}

/// The string form of a mapping key node, or `None` for non-scalar keys.
fn node_key(node: &MarkedYaml) -> Option<String> {
  match &node.data {
    YamlData::Value(Scalar::String(text)) => Some(text.to_string()),
    YamlData::Representation(text, _, _) => Some(text.to_string()),
    _ => None,
  }
}

/// Convert a node's char-offset span into a byte-offset miette span.
fn span_of(table: &[usize], node: &MarkedYaml) -> SourceSpan {
  let start = byte_of(table, node.span.start.index());
  let end = byte_of(table, node.span.end.index());
  SourceSpan::from((start, end.saturating_sub(start)))
}

/// A table mapping each char index to its byte offset, with a final sentinel at the text length.
fn char_byte_table(text: &str) -> Vec<usize> {
  let mut table: Vec<usize> = text.char_indices().map(|(byte, _)| byte).collect();
  table.push(text.len());
  table
}

/// Byte offset for a char index, clamped to the end of the text.
fn byte_of(table: &[usize], char_index: usize) -> usize {
  table
    .get(char_index)
    .copied()
    .unwrap_or_else(|| table.last().copied().unwrap_or(0))
}

/// Escape a mapping key for use as a JSON-pointer segment (RFC 6901).
pub(crate) fn escape_pointer(segment: &str) -> String {
  segment.replace('~', "~0").replace('/', "~1")
}

#[cfg(test)]
mod tests {
  use super::{ConfigSource, parse_config_document};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};

  fn parse(text: &str) -> Result<Value, Report> {
    let source = ConfigSource::new("config.yaml", text.to_owned());
    parse_config_document(&source, text)
  }

  fn parse_error_headline(text: &str) -> String {
    let err = parse(text).expect_err("expected a parse error");
    err.to_string().lines().next().unwrap_or_default().to_owned()
  }

  // A duplicate mapping key is a hard parse error (serde-saphyr `DuplicateKeyPolicy::Error`).
  #[test]
  fn test_source_parse_rejects_duplicate_mapping_key() {
    assert_eq!(
      "invalid configuration: could not parse config: error: line 2 column 1: duplicate mapping key: a, set DuplicateKeyPolicy in Options if acceptable",
      parse_error_headline("a: 1\na: 2\n")
    );
  }

  // A positive infinity literal in a typeless position is rejected (`reject_non_finite_typeless_float`).
  #[test]
  fn test_source_parse_rejects_infinity() {
    assert_eq!(
      "invalid configuration: could not parse config: error: line 1 column 4: non-finite float `.inf` rejected by reject_non_finite_typeless_float",
      parse_error_headline("x: .inf\n")
    );
  }

  // A not-a-number literal in a typeless position is rejected.
  #[test]
  fn test_source_parse_rejects_nan() {
    assert_eq!(
      "invalid configuration: could not parse config: error: line 1 column 4: non-finite float `.nan` rejected by reject_non_finite_typeless_float",
      parse_error_headline("x: .nan\n")
    );
  }

  // YAML 1.1 boolean words `no` and `on` deserialize as booleans, not strings.
  #[test]
  fn test_source_parse_yaml11_booleans_no_and_on() {
    let value = parse("first: no\nsecond: on\n").unwrap();
    assert_eq!(json!({ "first": false, "second": true }), value);
  }

  // A merge key (`<<`) is expanded into the surrounding mapping.
  #[test]
  fn test_source_parse_applies_merge_key() {
    let value = parse("base: &anchor\n  shared: 1\nchild:\n  <<: *anchor\n  own: 2\n").unwrap();
    assert_eq!(json!({ "shared": 1 }), value["base"]);
    assert_eq!(json!({ "shared": 1, "own": 2 }), value["child"]);
  }

  // A scientific-notation float keeps full f64 precision through the parse.
  #[test]
  fn test_source_parse_preserves_scientific_notation() {
    let value = parse("rate: 5.7e-05\n").unwrap();
    assert_eq!(json!({ "rate": 5.7e-05 }), value);
  }
}
