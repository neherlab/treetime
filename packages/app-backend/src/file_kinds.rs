use crate::command::{AppCommandExt, offered_commands};
use crate::runs::package_reader::PackageReader;
use app_commands::command::AppCommand;
use app_commands::config::properties::setting_properties;
use app_commands::yaml::yaml_value_read_str;
use app_datasets::text_schema_command;
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use strum_macros::Display;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime_io::csv::{DELIMITED_EXTENSIONS, default_metadata_delimiters, default_name_candidates};
use treetime_io::dates_csv::metadata_read_file;
use treetime_io::fasta::{FASTA_EXTENSIONS, fasta_read_file};
use treetime_io::nwk::TREE_DIALECT_DEFAULT;
use treetime_io::tree::{TREE_EXTENSIONS, tree_read_file};
use treetime_schema::skip_serializing_optionals;
use treetime_utils::error::report_to_string;
use treetime_utils::io::compression::remove_compression_ext;
use treetime_utils::io::file::read_file_with;
use treetime_utils::io::fs::extension;
use treetime_utils::{make_error, make_report};

const ZIP_SIGNATURES: [[u8; 4]; 2] = [*b"PK\x03\x04", *b"PK\x05\x06"];

const CONFIGURATION_EXTENSIONS: [&str; 3] = ["yaml", "yml", "json"];

const RUN_PACKAGE_EXTENSIONS: [&str; 1] = ["zip"];

const MIN_TREE_TIPS: usize = 2;

pub fn file_kinds(request: &FileKindsRequest) -> FileKinds {
  FileKinds {
    files: request.files.iter().map(|path| file_kind(path)).collect(),
  }
}

/// Files to identify by their content, for example files that the user dropped without naming the input they fill.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct FileKindsRequest {
  /// Paths of the files.
  pub files: Vec<PathBuf>,
}

/// The kind of each file of a request, in the order of the request.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct FileKinds {
  /// One entry per file of the request.
  pub files: Vec<DetectedFile>,
}

/// The kind of one file, as the readers of the commands read it.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct DetectedFile {
  /// Path of the file, as the request gives it.
  pub path: PathBuf,
  /// Kind of the file, from its content; the file name extension only decides which reader tries first.
  pub kind: FileKind,
  /// The analyses whose settings the configuration sets, for a `configuration`: the analysis that its
  /// `yaml-language-server` schema directive names, or else every analysis that has all of its top-level keys.
  pub commands: Vec<AppCommand>,
  /// Why no reader reads the file, for an `unknown` file.
  pub reason: Option<String>,
}

/// Kind of a file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, JsonSchema, Display, Serialize, Deserialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum FileKind {
  /// A tree in Newick or Nexus format, with at least two tips.
  Tree,
  /// Sequences in FASTA format.
  Alignment,
  /// A table with a column of sample names.
  Metadata,
  /// A YAML or JSON configuration of an analysis.
  Configuration,
  /// A zip archive whose `run.json` is the record of a run.
  RunPackage,
  /// None of the above.
  Unknown,
}

impl FileKind {
  const fn noun(self) -> &'static str {
    match self {
      Self::Tree => "a tree",
      Self::Alignment => "an alignment",
      Self::Metadata => "a metadata table",
      Self::Configuration => "a configuration",
      Self::RunPackage => "a run package",
      Self::Unknown => "a file of unknown kind",
    }
  }

  fn extensions(self) -> Vec<&'static str> {
    match self {
      Self::Tree => TREE_EXTENSIONS.to_vec(),
      Self::Alignment => FASTA_EXTENSIONS.to_vec(),
      Self::Metadata => DELIMITED_EXTENSIONS.iter().map(|(extension, _)| *extension).collect(),
      Self::Configuration => CONFIGURATION_EXTENSIONS.to_vec(),
      Self::RunPackage => RUN_PACKAGE_EXTENSIONS.to_vec(),
      Self::Unknown => vec![],
    }
  }
}

fn file_kind(path: &Path) -> DetectedFile {
  let unknown = |reason: String| DetectedFile {
    path: path.to_path_buf(),
    kind: FileKind::Unknown,
    commands: vec![],
    reason: Some(reason),
  };
  match is_zip(path) {
    Ok(true) => {
      return match read_run_package(path) {
        Ok(()) => DetectedFile {
          path: path.to_path_buf(),
          kind: FileKind::RunPackage,
          commands: vec![],
          reason: None,
        },
        Err(report) => unknown(sentence(&report_to_string(&report))),
      };
    },
    Ok(false) => {},
    Err(report) => return unknown(report_to_string(&report)),
  }
  let hint = extension_kind(path);
  let mut failures = vec![];
  for kind in reading_order(hint) {
    match read_as(kind, path) {
      Ok(commands) => {
        return DetectedFile {
          path: path.to_path_buf(),
          kind,
          commands,
          reason: None,
        };
      },
      Err(report) => failures.push(format!("as {}: {}", kind.noun(), report.root_cause())),
    }
  }
  let named = hint.map_or_else(String::new, |kind| {
    format!("The file name suggests {}, but ", kind.noun())
  });
  let reads = if named.is_empty() { "No reader" } else { "no reader" };
  unknown(format!(
    "{named}{reads} of TreeTime reads the file. Read {}",
    failures.join("; ")
  ))
}

fn reading_order(hint: Option<FileKind>) -> Vec<FileKind> {
  let order = [
    FileKind::Configuration,
    FileKind::Alignment,
    FileKind::Tree,
    FileKind::Metadata,
  ];
  hint
    .into_iter()
    .filter(|kind| order.contains(kind))
    .chain(order)
    .unique()
    .collect()
}

fn extension_kind(path: &Path) -> Option<FileKind> {
  let extension = extension(remove_compression_ext(path))?.to_lowercase();
  [
    FileKind::Tree,
    FileKind::Alignment,
    FileKind::Metadata,
    FileKind::Configuration,
    FileKind::RunPackage,
  ]
  .into_iter()
  .find(|kind| kind.extensions().contains(&extension.as_str()))
}

fn read_as(kind: FileKind, path: &Path) -> Result<Vec<AppCommand>, Report> {
  match kind {
    FileKind::Tree => read_tree(path).map(|()| vec![]),
    FileKind::Alignment => read_alignment(path).map(|()| vec![]),
    FileKind::Metadata => read_metadata(path).map(|()| vec![]),
    FileKind::Configuration => configuration_commands(path),
    FileKind::RunPackage | FileKind::Unknown => make_error!("TreeTime has no reader for {}", kind.noun()),
  }
}

fn read_tree(path: &Path) -> Result<(), Report> {
  let tree = tree_read_file(path, TREE_DIALECT_DEFAULT)?;
  let tips = tree.graph.num_leaves();
  if tips < MIN_TREE_TIPS {
    return make_error!("the tree has {tips} tips, and a tree for an analysis has at least {MIN_TREE_TIPS}");
  }
  Ok(())
}

fn read_alignment(path: &Path) -> Result<(), Report> {
  let records = match fasta_read_file(path, &Alphabet::new(AlphabetName::Nuc)?) {
    Ok(records) => records,
    Err(nucleotides) => match fasta_read_file(path, &Alphabet::new(AlphabetName::Aa)?) {
      Ok(records) => records,
      Err(amino_acids) => {
        let (nucleotides, amino_acids) = (
          nucleotides.root_cause().to_string(),
          amino_acids.root_cause().to_string(),
        );
        return if nucleotides == amino_acids {
          make_error!("{nucleotides}")
        } else {
          make_error!("as nucleotides, {nucleotides}; as amino acids, {amino_acids}")
        };
      },
    },
  };
  if records.is_empty() {
    return make_error!("the file holds no FASTA record");
  }
  Ok(())
}

fn read_metadata(path: &Path) -> Result<(), Report> {
  metadata_read_file(
    path,
    &default_metadata_delimiters(),
    &default_name_candidates(),
    None,
    None,
  )
  .map(|_| ())
}

fn configuration_commands(path: &Path) -> Result<Vec<AppCommand>, Report> {
  let text = read_file_with(path, |mut reader| {
    let mut text = String::new();
    reader
      .read_to_string(&mut text)
      .wrap_err("When reading the file as text")?;
    Ok(text)
  })?;
  if let Some(name) = text_schema_command(&text) {
    let command = AppCommand::iter()
      .find(|command| <&str>::from(*command) == name)
      .ok_or_else(|| make_report!("the schema directive names `{name}`, which is not a command of TreeTime"))?;
    if !command.in_apps() {
      return make_error!("the configuration is for `{command}`, which only the command-line tool runs");
    }
    return Ok(vec![command]);
  }
  let Value::Object(settings) = yaml_value_read_str(&text)? else {
    return make_error!("the document is not a mapping of settings");
  };
  let keys: BTreeSet<&str> = settings.keys().map(String::as_str).collect();
  let commands: Vec<AppCommand> = offered_commands()
    .map(|command| {
      let properties = setting_properties(command.config_schema().as_value())?;
      let known: BTreeSet<&str> = properties.iter().map(|property| property.key.as_str()).collect();
      Ok::<_, Report>(keys.is_subset(&known).then_some(command))
    })
    .filter_map_ok(|command| command)
    .try_collect()?;
  if commands.is_empty() {
    return make_error!(
      "no analysis has all of its settings: {}",
      keys.iter().map(|key| format!("`{key}`")).join(", ")
    );
  }
  Ok(commands)
}

fn is_zip(path: &Path) -> Result<bool, Report> {
  let mut file = File::open(path).wrap_err_with(|| format!("When opening file '{}'", path.display()))?;
  let mut signature = [0_u8; 4];
  let read = file
    .read(&mut signature)
    .wrap_err_with(|| format!("When reading file '{}'", path.display()))?;
  Ok(read == signature.len() && ZIP_SIGNATURES.contains(&signature))
}

fn read_run_package(path: &Path) -> Result<(), Report> {
  let file = File::open(path).wrap_err_with(|| format!("When opening file '{}'", path.display()))?;
  PackageReader::open(file)?.record().map(|_| ())
}

fn sentence(text: &str) -> String {
  let mut chars = text.chars();
  chars
    .next()
    .map_or_else(String::new, |first| first.to_uppercase().chain(chars).collect())
}
