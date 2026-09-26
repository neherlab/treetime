import type {
  AncestralConfig,
  AppCommand,
  ClockConfig,
  MugrationConfig,
  OptimizeConfig,
  PruneConfig,
  TimetreeConfig,
} from "@neherlab/app-contracts";

interface CommandConfigs {
  timetree: TimetreeConfig;
  clock: ClockConfig;
  ancestral: AncestralConfig;
  mugration: MugrationConfig;
  optimize: OptimizeConfig;
  prune: PruneConfig;
}

export type SettingKey<C extends AppCommand> = keyof CommandConfigs[C] & string;

interface CommandInfo {
  command: AppCommand;
  label: string;
  verb: string;
  description: string;
}

export const APP_COMMANDS: readonly AppCommand[] = ["timetree", "clock", "ancestral", "mugration", "optimize", "prune"];

export const COMMAND_INFO: Record<AppCommand, CommandInfo> = {
  timetree: {
    command: "timetree",
    label: "Time tree",
    verb: "Run time tree",
    description: "Date the ancestors and estimate the clock rate from sampling dates.",
  },
  clock: {
    command: "clock",
    label: "Clock signal",
    verb: "Run clock check",
    description: "Root-to-tip regression: clock rate, temporal signal and outliers.",
  },
  ancestral: {
    command: "ancestral",
    label: "Ancestral sequences",
    verb: "Reconstruct sequences",
    description: "Infer ancestral sequences and the mutations on each branch.",
  },
  mugration: {
    command: "mugration",
    label: "Discrete traits",
    verb: "Reconstruct traits",
    description: "Ancestral states of a metadata column, such as country or host.",
  },
  optimize: {
    command: "optimize",
    label: "Branch lengths",
    verb: "Optimize branch lengths",
    description: "Maximum-likelihood branch lengths on a fixed topology.",
  },
  prune: {
    command: "prune",
    label: "Prune tree",
    verb: "Prune tree",
    description: "Collapse short or empty branches and remove listed samples.",
  },
};

export const MAIN_SETTING_KEYS = {
  timetree: [
    "clock_rate",
    "clock_std_dev",
    "confidence",
    "covariation",
    "coalescent",
    "coalescent_opt",
    "coalescent_skyline",
    "skyline_n_points",
    "skyline_stiffness",
    "reroot",
    "keep_root",
    "clock_filter",
    "relax",
    "keep_polytomies",
    "model",
    "model_params",
    "max_iter",
  ],
  clock: ["reroot", "keep_root", "clock_filter", "covariation", "allow_negative_rate", "metadata_id_columns"],
  ancestral: ["method_anc", "model", "model_params", "gap_fill", "reconstruct_tip_states"],
  mugration: ["attribute", "weights", "pc", "missing_data", "sampling_bias_correction"],
  optimize: ["opt_method", "reroot", "divergence_units", "no_indels", "max_iter"],
  prune: ["prune_short", "prune_empty", "merge_shared_mutations", "prune_nodes_list"],
} satisfies { readonly [C in AppCommand]: ReadonlyArray<SettingKey<C>> };
