use serde::Deserialize;
use smart_default::SmartDefault;
use std::path::PathBuf;
use treetime::alphabet::alphabet::AlphabetName;
use treetime::ancestral::params::MethodAncestral;
use treetime::ancestral::sample::SampleMode;
use treetime::clock::find_best_root::params::RerootMethod;
use treetime::commands::ancestral::aa_model::AaModelName;
use treetime::commands::shared::reroot::RerootArgs;
use treetime::commands::timetree::args::TimeMarginalMode;
use treetime::gtr::get_gtr::GtrModelName;
use treetime::optimize::params::{BranchLengthMode, BranchOptMethod, InitialGuessMode};
use treetime::seq::gap_fill::GapFill;

use app_api::{
  TreetimeAncestralArgs, TreetimeClockArgs, TreetimeMugrationArgs, TreetimeOptimizeArgs, TreetimePruneArgs,
  TreetimeTimetreeArgs,
};

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct ServerAncestralArgs {
  pub input_fastas: Vec<String>,
  pub aln: Option<String>,
  pub vcf_reference: Option<String>,
  pub tree: String,
  pub alphabet: Option<AlphabetName>,
  #[default(GtrModelName::Infer)]
  pub model_name: GtrModelName,
  pub gtr_params: Vec<String>,
  #[default(MethodAncestral::default())]
  pub method_anc: MethodAncestral,
  pub dense: Option<bool>,
  pub aa: bool,
  #[default(GapFill::default())]
  pub gap_fill: GapFill,
  pub keep_overhangs: bool,
  pub zero_based: bool,
  pub reconstruct_tip_states: bool,
  pub report_ambiguous: bool,
  pub outdir: String,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: Option<u64>,
}

impl From<ServerAncestralArgs> for TreetimeAncestralArgs {
  fn from(s: ServerAncestralArgs) -> Self {
    use treetime::commands::shared::alignment::AlignmentArgs;
    use treetime::commands::shared::alphabet::AlphabetArgs;
    use treetime::commands::shared::gap_fill::GapFillArgs;
    use treetime::commands::shared::model::ModelArgs;
    use treetime::commands::shared::output::OutputArgs;
    Self {
      alignment: AlignmentArgs {
        alignment: s.input_fastas.into_iter().map(PathBuf::from).collect(),
      },
      vcf_reference: s.vcf_reference.map(PathBuf::from),
      tree: PathBuf::from(s.tree),
      alphabet_args: AlphabetArgs { alphabet: s.alphabet },
      model_args: ModelArgs {
        model: s.model_name,
        model_params: s.gtr_params,
      },
      method_anc: s.method_anc,
      dense: s.dense,
      gap_fill_args: GapFillArgs {
        gap_fill: s.gap_fill,
        keep_overhangs: s.keep_overhangs,
      },
      zero_based: s.zero_based,
      reconstruct_tip_states: s.reconstruct_tip_states,
      report_ambiguous: s.report_ambiguous,
      ignore_missing_alns: false,
      output_augur_node_data: None,
      translations: None,
      cdses: vec![],
      annotation: None,
      aa_root_sequence: None,
      aa_model: AaModelName::default(),
      output_aa_sequences: None,
      output: OutputArgs {
        outdir: PathBuf::from(s.outdir),
        ..Default::default()
      },
      gtr_iterations: s.gtr_iterations,
      site_specific_gtr: s.site_specific_gtr,
      seed: s.seed,
      aa: s.aa,
      marginal: false,
      custom_gtr: None,
      sample_from_profile: SampleMode::default(),
    }
  }
}

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct ServerClockArgs {
  pub aln: Vec<String>,
  pub tree: Option<String>,
  pub vcf_reference: Option<String>,
  pub dates: String,
  pub name_column: Option<String>,
  pub date_column: Option<String>,
  pub sequence_length: Option<usize>,
  #[default(GtrModelName::default())]
  pub gtr: GtrModelName,
  pub gtr_params: Vec<String>,
  #[default(BranchLengthMode::default())]
  pub branch_length_mode: BranchLengthMode,
  #[default(MethodAncestral::default())]
  pub method_anc: MethodAncestral,
  #[default = 3.0]
  pub clock_filter: f64,
  pub reroot: Option<RerootMethod>,
  pub reroot_tips: Vec<String>,
  pub keep_root: bool,
  pub prune_short: bool,
  pub tip_slack: Option<f64>,
  pub covariation: bool,
  pub allow_negative_rate: bool,
  pub outdir: String,
  pub seed: Option<u64>,
}

impl From<ServerClockArgs> for TreetimeClockArgs {
  fn from(s: ServerClockArgs) -> Self {
    use treetime::commands::shared::alignment::AlignmentArgs;
    use treetime::commands::shared::metadata::DateColumnArgs;
    use treetime::commands::shared::model::ModelArgs;
    use treetime::commands::shared::output::OutputArgs;
    Self {
      alignment: AlignmentArgs {
        alignment: s.aln.into_iter().map(PathBuf::from).collect(),
      },
      tree: s.tree.map(PathBuf::from),
      vcf_reference: s.vcf_reference.map(PathBuf::from),
      metadata: PathBuf::from(s.dates),
      metadata_id: {
        let mut id = treetime::commands::shared::metadata::MetadataIdArgs::default();
        if let Some(name_col) = s.name_column {
          id.metadata_id_columns = vec![name_col];
        }
        id
      },
      date_column: DateColumnArgs {
        date_column: s.date_column,
        ..DateColumnArgs::default()
      },
      sequence_length: s.sequence_length,
      model_args: ModelArgs {
        model: s.gtr,
        model_params: s.gtr_params,
      },
      branch_length_mode: s.branch_length_mode,
      method_anc: s.method_anc,
      clock_filter: s.clock_filter,
      reroot: RerootArgs {
        reroot: s.reroot,
        reroot_tips: s.reroot_tips,
      },
      keep_root: s.keep_root,
      prune_short: s.prune_short,
      tip_slack: s.tip_slack,
      covariation: s.covariation,
      allow_negative_rate: s.allow_negative_rate,
      output: OutputArgs {
        outdir: PathBuf::from(s.outdir),
        ..Default::default()
      },
      seed: s.seed,
      ..TreetimeClockArgs::default()
    }
  }
}

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct ServerTimetreeArgs {
  pub input_fastas: Vec<String>,
  pub tree: Option<String>,
  pub vcf_reference: Option<String>,
  pub dates: Option<String>,
  pub name_column: Option<String>,
  pub date_column: Option<String>,
  pub sequence_length: Option<usize>,
  pub clock_rate: Option<f64>,
  pub clock_std_dev: Option<f64>,
  #[default(BranchLengthMode::default())]
  pub branch_length_mode: BranchLengthMode,
  #[default(TimeMarginalMode::default())]
  pub time_marginal: TimeMarginalMode,
  pub confidence: bool,
  pub keep_polytomies: bool,
  pub resolve_polytomies: bool,
  pub relax: Vec<f64>,
  #[default = 2]
  pub max_iter: usize,
  pub coalescent: Option<f64>,
  pub coalescent_opt: bool,
  pub coalescent_skyline: bool,
  #[default = 10]
  pub n_skyline: usize,
  pub n_branches_posterior: Option<usize>,
  pub tip_labels: bool,
  pub no_tip_labels: bool,
  pub clock_filter: f64,
  pub n_iqd: Option<f64>,
  pub reroot: Option<RerootMethod>,
  pub reroot_tips: Vec<String>,
  pub keep_root: bool,
  pub allow_negative_rate: bool,
  pub tip_slack: Option<f64>,
  pub covariation: bool,
  #[default(GtrModelName::default())]
  pub gtr: GtrModelName,
  pub gtr_params: Vec<String>,
  #[default(MethodAncestral::default())]
  pub method_anc: MethodAncestral,
  #[default(AlphabetName::default())]
  pub alphabet: AlphabetName,
  pub dense: Option<bool>,
  pub aa: bool,
  #[default(GapFill::default())]
  pub gap_fill: GapFill,
  pub keep_overhangs: bool,
  pub zero_based: bool,
  pub reconstruct_tip_states: bool,
  pub report_ambiguous: bool,
  pub no_indels: bool,
  pub outdir: String,
  pub tracelog: Option<String>,
  pub seed: Option<u64>,
}

impl From<ServerTimetreeArgs> for TreetimeTimetreeArgs {
  fn from(s: ServerTimetreeArgs) -> Self {
    use treetime::commands::shared::alignment::AlignmentArgs;
    use treetime::commands::shared::alphabet::AlphabetArgs;
    use treetime::commands::shared::gap_fill::GapFillArgs;
    use treetime::commands::shared::metadata::DateColumnArgs;
    use treetime::commands::shared::model::ModelArgs;
    use treetime::commands::shared::output::{DivergenceUnits, OutputArgs};
    Self {
      alignment: AlignmentArgs {
        alignment: s.input_fastas.into_iter().map(PathBuf::from).collect(),
      },
      tree: s.tree.map(PathBuf::from),
      vcf_reference: s.vcf_reference.map(PathBuf::from),
      metadata: s.dates.map(PathBuf::from),
      metadata_id: {
        let mut id = treetime::commands::shared::metadata::MetadataIdArgs::default();
        if let Some(name_col) = s.name_column {
          id.metadata_id_columns = vec![name_col];
        }
        id
      },
      date_column_args: DateColumnArgs {
        date_column: s.date_column,
        ..DateColumnArgs::default()
      },
      sequence_length: s.sequence_length,
      clock_rate: s.clock_rate,
      clock_std_dev: s.clock_std_dev,
      branch_length_mode: s.branch_length_mode,
      time_marginal: s.time_marginal,
      confidence: s.confidence,
      keep_polytomies: s.keep_polytomies,
      resolve_polytomies: s.resolve_polytomies,
      relax: s.relax,
      max_iter: s.max_iter,
      coalescent: s.coalescent,
      coalescent_opt: s.coalescent_opt,
      coalescent_skyline: s.coalescent_skyline,
      n_skyline: s.n_skyline,
      n_branches_posterior: s.n_branches_posterior,
      plot_tree: None,
      plot_rtt: None,
      tip_labels: s.tip_labels,
      no_tip_labels: s.no_tip_labels,
      clock_filter: s.clock_filter,
      n_iqd: s.n_iqd,
      reroot: RerootArgs {
        reroot: s.reroot,
        reroot_tips: s.reroot_tips,
      },
      keep_root: s.keep_root,
      allow_negative_rate: s.allow_negative_rate,
      tip_slack: s.tip_slack,
      covariation: s.covariation,
      model_args: ModelArgs {
        model: s.gtr,
        model_params: s.gtr_params,
      },
      method_anc: s.method_anc,
      alphabet_args: AlphabetArgs {
        alphabet: Some(s.alphabet),
      },
      dense: s.dense,
      gap_fill_args: GapFillArgs {
        gap_fill: s.gap_fill,
        keep_overhangs: s.keep_overhangs,
      },
      zero_based: s.zero_based,
      reconstruct_tip_states: s.reconstruct_tip_states,
      report_ambiguous: s.report_ambiguous,
      no_indels: s.no_indels,
      divergence_units: DivergenceUnits::default(),
      output_augur_node_data: None,
      output: OutputArgs {
        outdir: PathBuf::from(s.outdir),
        ..Default::default()
      },
      tracelog: s.tracelog.map(PathBuf::from),
      seed: s.seed,
      aa: s.aa,
      custom_gtr: None,
      clock_filter_method: None,
      gen_per_year: None,
      greedy_resolve: false,
      stochastic_resolve: false,
    }
  }
}

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct ServerMugrationArgs {
  pub tree: Option<String>,
  #[default(_code = r#""country".to_owned()"#)]
  pub attribute: String,
  pub states: String,
  pub weights: Option<String>,
  pub name_column: Option<String>,
  pub confidence: Option<String>,
  pub pc: Option<f64>,
  #[default(_code = r#""?".to_owned()"#)]
  pub missing_data: String,
  #[default = 0.5]
  pub missing_weights_threshold: f64,
  #[default = 5]
  pub iterations: usize,
  pub sampling_bias_correction: Option<f64>,
  pub outdir: String,
}

impl From<ServerMugrationArgs> for TreetimeMugrationArgs {
  fn from(s: ServerMugrationArgs) -> Self {
    use treetime::commands::shared::output::OutputArgs;
    Self {
      tree: s.tree.map(PathBuf::from),
      attribute: s.attribute,
      metadata: PathBuf::from(s.states),
      weights: s.weights.map(PathBuf::from),
      metadata_id: {
        let mut id = treetime::commands::shared::metadata::MetadataIdArgs::default();
        if let Some(name_col) = s.name_column {
          id.metadata_id_columns = vec![name_col];
        }
        id
      },
      output_confidence: s.confidence.map(PathBuf::from),
      pc: s.pc,
      missing_data: s.missing_data,
      missing_weights_threshold: s.missing_weights_threshold,
      iterations: s.iterations,
      sampling_bias_correction: s.sampling_bias_correction,
      smooth_initial_pi: false,
      filter_uninformative_root: false,
      output_augur_node_data: None,
      seed: None,
      output: OutputArgs {
        outdir: PathBuf::from(s.outdir),
        ..Default::default()
      },
    }
  }
}

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct ServerOptimizeArgs {
  pub input_fastas: Vec<String>,
  pub tree: String,
  pub alphabet: Option<AlphabetName>,
  #[default(GtrModelName::Infer)]
  pub model_name: GtrModelName,
  pub dense: Option<bool>,
  pub outdir: String,
  #[default = 10]
  pub max_iter: usize,
  #[default = 0.1]
  pub dp: f64,
  #[default = 0.75]
  pub damping: f64,
  #[default(InitialGuessMode::Auto)]
  pub branch_length_initial_guess: InitialGuessMode,
  #[default(BranchOptMethod::default())]
  pub opt_method: BranchOptMethod,
  pub no_indels: bool,
  #[default(GapFill::default())]
  pub gap_fill: GapFill,
  pub keep_overhangs: bool,
}

impl From<ServerOptimizeArgs> for TreetimeOptimizeArgs {
  fn from(s: ServerOptimizeArgs) -> Self {
    use treetime::commands::shared::alignment::AlignmentArgs;
    use treetime::commands::shared::alphabet::AlphabetArgs;
    use treetime::commands::shared::gap_fill::GapFillArgs;
    use treetime::commands::shared::model::ModelArgs;
    use treetime::commands::shared::output::{DivergenceUnits, OutputArgs};
    Self {
      alignment: AlignmentArgs {
        alignment: s.input_fastas.into_iter().map(PathBuf::from).collect(),
      },
      tree: PathBuf::from(s.tree),
      alphabet_args: AlphabetArgs { alphabet: s.alphabet },
      model_args: ModelArgs {
        model: s.model_name,
        ..ModelArgs::default()
      },
      dense: s.dense,
      output: OutputArgs {
        outdir: PathBuf::from(s.outdir),
        ..Default::default()
      },
      max_iter: s.max_iter,
      dp: s.dp,
      damping: s.damping,
      branch_length_initial_guess: s.branch_length_initial_guess,
      opt_method: s.opt_method,
      no_indels: s.no_indels,
      divergence_units: DivergenceUnits::default(),
      gap_fill_args: GapFillArgs {
        gap_fill: s.gap_fill,
        keep_overhangs: s.keep_overhangs,
      },
      output_augur_node_data: None,
    }
  }
}

#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct ServerPruneArgs {
  pub input_fastas: Vec<String>,
  pub tree: String,
  pub alphabet: Option<AlphabetName>,
  pub outdir: String,
  pub prune_short: Option<f64>,
  pub prune_empty: bool,
  pub merge_shared_mutations: bool,
  pub prune_nodes_list: Option<String>,
  #[default = ',']
  pub prune_nodes_list_delimiter: char,
  pub prune_nodes_list_file: Option<String>,
  #[default = '\n']
  pub prune_nodes_list_file_delimiter: char,
}

impl From<ServerPruneArgs> for TreetimePruneArgs {
  fn from(s: ServerPruneArgs) -> Self {
    use treetime::commands::shared::alignment::AlignmentArgs;
    use treetime::commands::shared::alphabet::AlphabetArgs;
    use treetime::commands::shared::output::OutputArgs;
    Self {
      alignment: AlignmentArgs {
        alignment: s.input_fastas.into_iter().map(PathBuf::from).collect(),
      },
      tree: PathBuf::from(s.tree),
      alphabet_args: AlphabetArgs { alphabet: s.alphabet },
      output: OutputArgs {
        outdir: PathBuf::from(s.outdir),
        ..Default::default()
      },
      prune_short: s.prune_short,
      prune_empty: s.prune_empty,
      merge_shared_mutations: s.merge_shared_mutations,
      prune_nodes_list: s.prune_nodes_list,
      prune_nodes_list_delimiter: s.prune_nodes_list_delimiter,
      prune_nodes_list_file: s.prune_nodes_list_file.map(PathBuf::from),
      prune_nodes_list_file_delimiter: s.prune_nodes_list_file_delimiter,
    }
  }
}
