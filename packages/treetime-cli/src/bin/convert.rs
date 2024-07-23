use clap::{AppSettings, ArgEnum, Parser, ValueHint};
use clap_verbosity_flag::{Verbosity, WarnLevel};
use color_eyre::{Section, SectionExt};
use ctor::ctor;
use eyre::Report;
use lazy_static::lazy_static;
use log::LevelFilter;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use treetime::graph::edge::{GraphEdge, Weighted};
use treetime::graph::graph::Graph;
use treetime::graph::node::{GraphNode, Named};
use treetime::io::auspice::{
  auspice_read_file, auspice_write_file, AuspiceDataFromGraphData, AuspiceDataToGraphData, AuspiceFromGraph,
  AuspiceGraphContext, AuspiceToGraph, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeContext, AuspiceTreeData,
  AuspiceTreeMeta, AuspiceTreeNode, AuspiceTreeNodeAttrs,
};
use treetime::io::compression::remove_compression_ext;
use treetime::io::fs::extension;
use treetime::io::graphviz::{EdgeToGraphViz, NodeToGraphviz};
use treetime::io::json::{json_read_file, json_write_file, JsonPretty};
use treetime::io::nex::{nex_write_file, NexWriteOptions};
use treetime::io::nwk::{
  format_weight, nwk_read_file, nwk_write_file, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions,
};
use treetime::io::usher_mat::{
  usher_mat_pb_read_file, usher_mat_pb_write_file, UsherDataToGraphData, UsherFromGraph, UsherGraphContext,
  UsherMetadata, UsherMutationList, UsherToGraph, UsherTree, UsherTreeContext, UsherTreeNode,
};
use treetime::utils::global_init::{global_init, setup_logger};
use treetime::{make_internal_report, make_report};

#[ctor]
fn init() {
  global_init();
}

fn main() -> Result<(), Report> {
  rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

  let args = Args::parse();

  // --verbosity=<level> and --silent take priority over -v and -q
  let filter_level = if args.silent {
    LevelFilter::Off
  } else {
    match args.verbosity {
      None => args.verbose.log_level_filter(),
      Some(verbosity) => verbosity,
    }
  };

  setup_logger(filter_level);

  let input_format = args
    .input_format
    .or_else(|| guess_tree_format_from_filename(remove_compression_ext(&args.input)))
    .ok_or_else(|| {
      make_report!("Input format was not specified and unable to autodetect. Please provide --input-format argument")
    })
    .with_section(|| format!("{:#?}", &args.input).header("Input file:"))?;

  let output_format = args
    .output_format
    .or_else(|| guess_tree_format_from_filename(remove_compression_ext(&args.output)))
    .ok_or_else(|| {
      make_report!("Output format was not specified and unable to autodetect. Please provide --output-format argument")
    })
    .with_section(|| format!("{:#?}", &args.input).header("Output file:"))?;

  let graph: ConverterGraph = match input_format {
    TreeFormat::Auspice => auspice_read_file(&args.input),
    TreeFormat::Newick => nwk_read_file(&args.input),
    TreeFormat::Nexus => unimplemented!("Reading Nexus files is not yet implemented"),
    TreeFormat::PhyloGraph => json_read_file(&args.input),
    // TreeFormat::MatJson => usher_mat_json_read_file(&args.input),
    TreeFormat::MatPb => usher_mat_pb_read_file(&args.input),
    _ => panic!(),
  }?;

  match output_format {
    TreeFormat::Auspice => auspice_write_file(&args.output, &graph),
    TreeFormat::Newick => nwk_write_file(&args.output, &graph, &NwkWriteOptions::default()),
    TreeFormat::Nexus => nex_write_file(&args.output, &graph, &NexWriteOptions::default()),
    TreeFormat::PhyloGraph => json_write_file(&args.output, &graph, JsonPretty(true)),
    // TreeFormat::MatJson => usher_mat_json_write_file(&args.output, &graph),
    TreeFormat::MatPb => usher_mat_pb_write_file(&args.output, &graph),
    _ => panic!(),
  }?;

  Ok(())
}

type ConverterGraph = Graph<ConverterNode, ConverterEdge, ConverterData>;

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ConverterNode {
  pub name: Option<String>,
}

impl GraphNode for ConverterNode {}

impl Named for ConverterNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl NodeFromNwk for ConverterNode {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|n| n.as_ref().to_owned()),
    })
  }
}

impl NodeToNwk for ConverterNode {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

impl NodeToGraphviz for ConverterNode {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ConverterEdge {
  pub weight: Option<f64>,
}

impl GraphEdge for ConverterEdge {}

impl Weighted for ConverterEdge {
  fn weight(&self) -> Option<f64> {
    self.weight
  }
  fn set_weight(&mut self, weight: Option<f64>) {
    self.weight = weight;
  }
}

impl EdgeFromNwk for ConverterEdge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self { weight })
  }
}

impl EdgeToNwk for ConverterEdge {
  fn nwk_weight(&self) -> Option<f64> {
    self.weight
  }
}

impl EdgeToGraphViz for ConverterEdge {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .weight
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.weight
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct ConverterData {
  pub version: Option<String>,
  pub meta: AuspiceTreeMeta,
  pub root_sequence: Option<BTreeMap<String, String>>,
  pub other: Value,
}

impl AuspiceDataFromGraphData for ConverterData {
  fn auspice_data_from_graph_data(&self) -> Result<AuspiceTreeData, Report> {
    Ok(AuspiceTreeData {
      version: self.version.clone(),
      meta: self.meta.clone(),
      root_sequence: self.root_sequence.clone(),
      other: self.other.clone(),
    })
  }
}

impl AuspiceDataToGraphData for ConverterData {
  fn auspice_data_to_graph_data(tree: &AuspiceTree) -> Result<Self, Report> {
    Ok(Self {
      version: tree.data.version.clone(),
      meta: tree.data.meta.clone(),
      root_sequence: tree.data.root_sequence.clone(),
      other: tree.data.other.clone(),
    })
  }
}

impl AuspiceFromGraph<ConverterNode, ConverterEdge, ConverterData> for () {
  fn auspice_node_from_graph_components(
    AuspiceGraphContext { node, edge, .. }: &AuspiceGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<AuspiceTreeNode, Report> {
    Ok(AuspiceTreeNode {
      name: node.name.clone().unwrap_or_default(),
      branch_attrs: AuspiceTreeBranchAttrs::default(),
      node_attrs: AuspiceTreeNodeAttrs {
        div: edge.and_then(|edge| edge.weight),
        clade_membership: None,
        region: None,
        country: None,
        division: None,
        other: Value::default(),
      },
      children: vec![],
      other: Value::default(),
    })
  }
}

impl AuspiceToGraph<ConverterNode, ConverterEdge, ConverterData> for () {
  fn auspice_node_to_graph_components(
    AuspiceTreeContext { node, .. }: &AuspiceTreeContext,
  ) -> Result<(ConverterNode, ConverterEdge), Report> {
    Ok((
      ConverterNode {
        name: Some(node.name.clone()),
      },
      ConverterEdge {
        weight: node.node_attrs.div,
      },
    ))
  }
}

impl UsherDataToGraphData for ConverterData {
  fn usher_data_to_graph_data(_: &UsherTree) -> Result<Self, Report> {
    Ok(Self::default())
  }
}

impl UsherFromGraph<ConverterNode, ConverterEdge, ConverterData> for () {
  fn usher_node_from_graph_components(
    context: &UsherGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report> {
    let &UsherGraphContext { node, .. } = context;

    let node_name = node
      .name
      .as_ref()
      .ok_or_else(|| make_internal_report!("Encountered node with empty name"))?
      .to_owned();

    let usher_node = UsherTreeNode {
      node_name,
      condensed_leaves: vec![],
    };

    let usher_mutations = UsherMutationList { mutation: vec![] };

    let usher_meta = UsherMetadata {
      clade_annotations: vec![],
    };

    Ok((usher_node, usher_mutations, usher_meta))
  }
}

impl UsherToGraph<ConverterNode, ConverterEdge, ConverterData> for () {
  fn usher_node_to_graph_components(context: &UsherTreeContext) -> Result<(ConverterNode, ConverterEdge), Report> {
    let UsherTreeContext { node, tree } = context;

    let node = ConverterNode {
      name: node.name.as_ref().map(ToOwned::to_owned),
    };

    let edge = ConverterEdge { weight: None };

    Ok((node, edge))
  }
}

#[derive(Copy, Clone, Debug, ArgEnum)]
#[clap(rename = "kebab-case")]
pub enum TreeFormat {
  Auspice,
  MatJson,
  MatPb,
  Newick,
  Nexus,
  PhyloGraph,
}

pub fn guess_tree_format_from_filename(filepath: impl AsRef<Path>) -> Option<TreeFormat> {
  let filepath = filepath.as_ref();
  let ext = extension(filepath).map(|s| s.to_lowercase());
  match ext.as_deref() {
    Some("auspice.json") => Some(TreeFormat::Auspice),
    Some("graph.json") => Some(TreeFormat::PhyloGraph),
    Some("mat.json") => Some(TreeFormat::MatPb),
    Some("mat.pb") => Some(TreeFormat::MatPb),
    Some("nex | nexus") => Some(TreeFormat::Nexus),
    Some("nwk" | "newick") => Some(TreeFormat::Newick),
    _ => None,
  }
}

#[derive(Parser, Debug)]
#[clap(name = "treetime", trailing_var_arg = true)]
#[clap(author, version)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(verbatim_doc_comment)]
/// Read and write Usher MAT files
///
/// * https://github.com/yatisht/usher
/// * https://usher-wiki.readthedocs.io/
pub struct Args {
  /// Path to input file
  ///
  /// Accepts plain or compressed files. If a compressed file is provided, it will be transparently
  /// decompressed. Supported compression formats: "gz", "bz2", "xz", "zst". Decompressor is chosen based on file
  /// extension.
  ///
  /// Omit this argument or use special value "-" to read uncompressed data from standard input (stdin).
  #[clap(value_hint = ValueHint::FilePath)]
  #[clap(display_order = 1)]
  #[clap(default_value = "-")]
  pub input: PathBuf,

  /// Path to output file
  ///
  /// If the provided file path ends with one of the supported extensions: "gz", "bz2", "xz", "zst", then the file will be written compressed. Compressor is chosen based on file extension.
  ///
  /// Omit this argument or use special value "-" to write uncompressed data to standard output (stdout).
  #[clap(long, short = 'o')]
  #[clap(display_order = 1)]
  #[clap(value_hint = ValueHint::AnyPath)]
  #[clap(default_value = "-")]
  pub output: PathBuf,

  /// Input format to read
  ///
  /// Provide this argument if automatic format detection fails or if you want to override it.
  #[clap(long, short = 'r', arg_enum)]
  #[clap(display_order = 2)]
  pub input_format: Option<TreeFormat>,

  /// Output format to write
  ///
  /// Provide this argument if automatic format detection fails or if you want to override it.
  #[clap(long, short = 'w', arg_enum)]
  #[clap(display_order = 2)]
  pub output_format: Option<TreeFormat>,

  /// Make output more quiet or more verbose
  #[clap(flatten)]
  pub verbose: Verbosity<WarnLevel>,

  /// Set verbosity level
  #[clap(long, global = true, conflicts_with = "verbose", conflicts_with = "silent", possible_values(VERBOSITIES.iter()))]
  #[clap(display_order = 999)]
  pub verbosity: Option<LevelFilter>,

  /// Disable all console output. Same as --verbosity=off
  #[clap(long, global = true, conflicts_with = "verbose", conflicts_with = "verbosity")]
  #[clap(display_order = 999)]
  pub silent: bool,
}

lazy_static! {
  static ref VERBOSITIES: &'static [&'static str] = &["off", "error", "warn", "info", "debug", "trace"];
}
