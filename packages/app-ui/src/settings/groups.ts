import type { SettingSpec } from "./schema";

const SETTING_GROUPS = [
  "Input data",
  "Molecular clock",
  "Rooting",
  "Dating",
  "Coalescent prior",
  "Polytomies",
  "Substitution model",
  "Ancestral reconstruction",
  "Branch lengths",
  "Tree order",
  "Outputs",
  "Reproducibility",
  "Other",
] as const;

export type SettingGroup = (typeof SETTING_GROUPS)[number];

const GROUP_KEYS: ReadonlyArray<readonly [SettingGroup, readonly string[]]> = [
  [
    "Input data",
    [
      "tree",
      "alignment",
      "metadata",
      "vcf_reference",
      "metadata_id_columns",
      "metadata_delimiters",
      "date_column",
      "date_format",
      "sequence_length",
      "attribute",
      "weights",
      "missing_data",
      "annotation",
      "cdses",
      "translations",
      "prune_nodes_list",
      "prune_nodes_list_delimiter",
      "prune_nodes_list_file",
      "prune_nodes_list_file_delimiter",
      "ignore_missing_alns",
      "zero_based",
      "alphabet",
      "aa",
    ],
  ],
  [
    "Molecular clock",
    [
      "clock_rate",
      "clock_std_dev",
      "covariation",
      "clock_filter",
      "n_iqd",
      "clock_filter_method",
      "prune_outliers",
      "allow_negative_rate",
      "tip_slack",
      "relax",
      "clock_regression",
    ],
  ],
  ["Rooting", ["reroot", "reroot_tips", "keep_root", "branch_split"]],
  [
    "Dating",
    ["confidence", "time_marginal", "branch_length_mode", "max_iter", "n_branches_posterior", "divergence_units"],
  ],
  [
    "Coalescent prior",
    [
      "coalescent",
      "coalescent_opt",
      "coalescent_skyline",
      "skyline_n_points",
      "skyline_stiffness",
      "coalescent_confidence",
      "gen_per_year",
    ],
  ],
  ["Polytomies", ["keep_polytomies", "resolve_polytomies", "greedy_resolve", "stochastic_resolve"]],
  [
    "Substitution model",
    [
      "model",
      "model_params",
      "custom_gtr",
      "aa_model",
      "site_specific_gtr",
      "gtr_iterations",
      "smooth_initial_pi",
      "pc",
      "iterations",
      "missing_weights_threshold",
      "sampling_bias_correction",
      "filter_uninformative_root",
    ],
  ],
  [
    "Ancestral reconstruction",
    [
      "method_anc",
      "marginal",
      "dense",
      "gap_fill",
      "keep_overhangs",
      "include_leaves",
      "impute_missing_data",
      "reconstruct_tip_states",
      "report_ambiguous",
      "no_indels",
      "sample_from_profile",
      "aa_root_sequence",
    ],
  ],
  [
    "Branch lengths",
    [
      "opt_method",
      "branch_length_initial_guess",
      "damping",
      "dp",
      "no_collapse_short_branches",
      "no_flip_parent_child",
      "no_merge_siblings",
      "prune_short",
      "prune_empty",
      "merge_shared_mutations",
    ],
  ],
  [
    "Tree order",
    [
      "ladderize",
      "topology_order",
      "topology_order_target_source",
      "topology_order_target_file",
      "topology_order_target_aggregate",
      "tip_labels",
      "no_tip_labels",
    ],
  ],
  ["Outputs", ["output_selection", "output_nwk_style", "plot_tree", "plot_rtt"]],
  ["Reproducibility", ["seed"]],
];

const GROUP_OF_KEY = new Map(GROUP_KEYS.flatMap(([group, keys]) => keys.map((key) => [key, group] as const)));

const ORDER_OF_KEY = new Map(GROUP_KEYS.flatMap(([, keys]) => keys.map((key, index) => [key, index] as const)));

export function groupedSpecs(specs: readonly SettingSpec[]): Array<[SettingGroup, SettingSpec[]]> {
  const sorted = specs.toSorted(
    (left, right) =>
      (ORDER_OF_KEY.get(left.path[0] ?? "") ?? Number.MAX_SAFE_INTEGER) -
        (ORDER_OF_KEY.get(right.path[0] ?? "") ?? Number.MAX_SAFE_INTEGER) || left.key.localeCompare(right.key),
  );

  return SETTING_GROUPS.flatMap((group) => {
    const members = sorted.filter((spec) => groupOf(spec) === group);

    return members.length > 0 ? [[group, members]] : [];
  });
}

function groupOf(spec: SettingSpec): SettingGroup {
  if (spec.pathRole === "output") {
    return "Outputs";
  }

  return GROUP_OF_KEY.get(spec.path[0] ?? "") ?? "Other";
}
