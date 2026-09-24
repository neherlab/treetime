use bon::bon;
use color_eyre::Section;
use eyre::Report;
use itertools::Itertools;
use miette::{Diagnostic, LabeledSpan, NamedSource, Severity, SourceCode, SourceSpan};
use saphyr::{LoadableYamlNode, MarkedYaml, Scalar, YamlData};
use serde_json::Value;
use serde_saphyr::DuplicateKeyPolicy;
use std::collections::BTreeMap;
use std::fmt::Display;
use treetime_utils::{make_error, make_report};

pub(crate) fn parse_config_document(source: &ConfigSource, text: &str) -> Result<Value, Report> {
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
        vec![RawDiagnostic::builder("config::syntax", format!("could not parse config: {err}")).build()],
      )?;
      unreachable!("render_and_bail returns an error whenever diagnostics are present");
    },
  }
}

pub(crate) fn render_and_bail(
  source: &ConfigSource,
  top_message: &str,
  diags: Vec<RawDiagnostic>,
) -> Result<(), Report> {
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

pub(crate) struct RawDiagnostic {
  pub pointer: Option<String>,
  pub use_key_span: bool,
  pub code: String,
  pub message: String,
  pub label: String,
  pub help: Option<String>,
}

#[bon]
impl RawDiagnostic {
  #[builder]
  pub(crate) fn new(
    #[builder(start_fn, into)] code: String,
    #[builder(start_fn, into)] message: String,
    #[builder(into, name = at)] pointer: Option<String>,
    #[builder(default, name = key_span)] use_key_span: bool,
    #[builder(into)] help: Option<String>,
  ) -> Self {
    Self {
      pointer,
      use_key_span,
      code,
      message,
      label: "here".to_owned(),
      help,
    }
  }

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

pub(crate) struct ConfigSource {
  name: String,
  text: String,
  spans: BTreeMap<String, NodeSpan>,
}

impl ConfigSource {
  #[cfg_attr(
    dylint_lib = "treetime_lints",
    expect(
      error_dropped_by_pattern,
      reason = "source spans are optional; parse_config_document reports the YAML error for the same text"
    )
  )]
  pub(crate) fn new(name: impl Into<String>, text: impl Into<String>) -> Self {
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

  fn span_for(&self, pointer: &str) -> Option<SourceSpan> {
    self.spans.get(pointer).map(|node| node.value)
  }

  fn key_span_for(&self, pointer: &str) -> Option<SourceSpan> {
    self.spans.get(pointer).map(|node| node.key.unwrap_or(node.value))
  }

  fn named_source(&self) -> NamedSource<String> {
    NamedSource::new(&self.name, self.text.clone())
  }
}

#[derive(Debug, derive_more::Display, derive_more::Error)]
#[display("{message}")]
struct ConfigReport {
  message: String,
  related: Vec<ConfigDiagnostic>,
}

impl Diagnostic for ConfigReport {
  fn severity(&self) -> Option<Severity> {
    Some(Severity::Error)
  }

  fn related<'a>(&'a self) -> Option<Box<dyn Iterator<Item = &'a dyn Diagnostic> + 'a>> {
    Some(Box::new(self.related.iter().map(|diag| -> &dyn Diagnostic { diag })))
  }
}

#[derive(Debug, derive_more::Display, derive_more::Error)]
#[display("{message}")]
struct ConfigDiagnostic {
  message: String,
  code: String,
  help: Option<String>,
  src: NamedSource<String>,
  span: Option<SourceSpan>,
  label: String,
}

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

struct NodeSpan {
  value: SourceSpan,
  key: Option<SourceSpan>,
}

fn node_key(node: &MarkedYaml) -> Option<String> {
  match &node.data {
    YamlData::Value(Scalar::String(text)) => Some(text.to_string()),
    YamlData::Representation(text, _, _) => Some(text.to_string()),
    _ => None,
  }
}

fn span_of(table: &[usize], node: &MarkedYaml) -> SourceSpan {
  let start = byte_of(table, node.span.start.index());
  let end = byte_of(table, node.span.end.index());
  SourceSpan::from((start, end.saturating_sub(start)))
}

fn char_byte_table(text: &str) -> Vec<usize> {
  let mut table: Vec<usize> = text.char_indices().map(|(byte, _)| byte).collect();
  table.push(text.len());
  table
}

fn byte_of(table: &[usize], char_index: usize) -> usize {
  table
    .get(char_index)
    .copied()
    .unwrap_or_else(|| table.last().copied().unwrap_or(0))
}

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

  #[test]
  fn test_source_parse_rejects_duplicate_mapping_key() {
    assert_eq!(
      "invalid configuration: could not parse config: error: line 2 column 1: duplicate mapping key: a, set DuplicateKeyPolicy in Options if acceptable",
      parse_error_headline("a: 1\na: 2\n")
    );
  }

  #[test]
  fn test_source_parse_rejects_infinity() {
    assert_eq!(
      "invalid configuration: could not parse config: error: line 1 column 4: non-finite float `.inf` rejected by reject_non_finite_typeless_float",
      parse_error_headline("x: .inf\n")
    );
  }

  #[test]
  fn test_source_parse_rejects_nan() {
    assert_eq!(
      "invalid configuration: could not parse config: error: line 1 column 4: non-finite float `.nan` rejected by reject_non_finite_typeless_float",
      parse_error_headline("x: .nan\n")
    );
  }

  #[test]
  fn test_source_parse_yaml11_booleans_no_and_on() {
    let value = parse("first: no\nsecond: on\n").unwrap();
    assert_eq!(json!({ "first": false, "second": true }), value);
  }

  #[test]
  fn test_source_parse_applies_merge_key() {
    let value = parse("base: &anchor\n  shared: 1\nchild:\n  <<: *anchor\n  own: 2\n").unwrap();
    assert_eq!(json!({ "shared": 1 }), value["base"]);
    assert_eq!(json!({ "shared": 1, "own": 2 }), value["child"]);
  }

  #[test]
  fn test_source_parse_preserves_scientific_notation() {
    let value = parse("rate: 5.7e-05\n").unwrap();
    assert_eq!(json!({ "rate": 5.7e-05 }), value);
  }
}
