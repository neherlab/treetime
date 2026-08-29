use eyre::Report;
use miette::{Diagnostic, LabeledSpan, NamedSource, Severity, SourceCode, SourceSpan};
use saphyr::{LoadableYamlNode, MarkedYaml, Scalar, YamlData};
use std::collections::BTreeMap;
use std::fmt::{self, Display, Formatter};
use treetime_utils::make_error;

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

  /// Replace the default caret label.
  #[must_use]
  pub fn label(mut self, label: impl Into<String>) -> Self {
    self.label = label.into();
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

/// Render all diagnostics as one batched miette report to stderr and return a terse eyre error.
///
/// The detailed, caret-annotated report is printed here; the returned error is a single line so that
/// the globally installed color-eyre hook does not print the details a second time (D13). An empty
/// diagnostic list is success.
pub fn render_and_bail(source: &ConfigSource, top_message: &str, diags: Vec<RawDiagnostic>) -> Result<(), Report> {
  if diags.is_empty() {
    return Ok(());
  }
  let count = diags.len();
  let related = diags.into_iter().map(|diag| diag.resolve(source)).collect();
  let report = ConfigReport {
    message: top_message.to_owned(),
    related,
  };

  let mut rendered = String::new();
  if miette::GraphicalReportHandler::new()
    .render_report(&mut rendered, &report)
    .is_ok()
  {
    eprint!("{rendered}");
  }

  let plural = if count == 1 { "" } else { "s" };
  make_error!("{top_message} ({count} problem{plural} reported above)")
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
fn escape_pointer(segment: &str) -> String {
  segment.replace('~', "~0").replace('/', "~1")
}
