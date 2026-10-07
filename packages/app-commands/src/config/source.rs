use crate::yaml::yaml_value_read_str;
use bon::bon;
use color_eyre::Section;
use deser::Serialize;
use eyre::Report;
use itertools::Itertools;
use miette::{Diagnostic, LabeledSpan, NamedSource, Severity, SourceCode, SourceSpan};
use saphyr::{LoadableYamlNode, MarkedYaml, Scalar, YamlData};
use schemars::JsonSchema;
use serde_json::Value;
use std::collections::BTreeMap;
use std::fmt::Display;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::make_report;

pub fn parse_config_document(source: &ConfigSource, text: &str) -> Result<Value, Report> {
  match yaml_value_read_str(text) {
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
  let rendered = rendered.trim_end().to_owned();

  let invalid = InvalidConfig {
    message: format!("{top_message}: {problems}"),
    problems: report.related.iter().map(ConfigDiagnostic::problem).collect(),
    rendered: rendered.clone(),
  };
  Err(Report::new(invalid).with_section(move || rendered))
}

/// Configuration rejected by parsing or by the schema check.
#[derive(Clone, Debug, derive_more::Display, derive_more::Error, JsonSchema)]
#[display("{message}")]
pub struct InvalidConfig {
  /// One-line summary of every problem, as the CLI prints it.
  pub message: String,
  /// Each problem with its location in the config.
  pub problems: Vec<ConfigProblem>,
  /// The problems drawn against the config text, as the CLI prints them.
  pub rendered: String,
}

impl InvalidConfig {
  pub fn problems_of(report: &Report) -> Vec<ConfigProblem> {
    report
      .downcast_ref::<Self>()
      .map(|invalid| invalid.problems.clone())
      .unwrap_or_default()
  }
}

/// One problem found in a configuration.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct ConfigProblem {
  /// Stable diagnostic code, for example `config::unknown-field`.
  pub code: String,
  /// Human-readable description of the problem.
  pub message: String,
  /// Byte offset and length of the offending text in the config, when known.
  pub span: Option<ConfigSpan>,
  /// Suggestion for fixing the problem.
  pub help: Option<String>,
}

/// Location of a problem in the configuration text.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct ConfigSpan {
  /// Byte offset of the first character.
  pub offset: usize,
  /// Length in bytes.
  pub length: usize,
}

pub struct RawDiagnostic {
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
  pub fn new(
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

pub struct ConfigSource {
  name: String,
  text: String,
  spans: BTreeMap<String, NodeSpan>,
}

impl ConfigSource {
  #[cfg_attr(
    dylint_lib = "custom",
    expect(
      error_dropped_by_pattern,
      reason = "source spans are optional; parse_config_document reports the YAML error for the same text"
    )
  )]
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

impl ConfigDiagnostic {
  fn problem(&self) -> ConfigProblem {
    ConfigProblem {
      code: self.code.clone(),
      message: self.message.clone(),
      span: self.span.map(|span| ConfigSpan {
        offset: span.offset(),
        length: span.len(),
      }),
      help: self.help.clone(),
    }
  }
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

pub fn escape_pointer(segment: &str) -> String {
  segment.replace('~', "~0").replace('/', "~1")
}
