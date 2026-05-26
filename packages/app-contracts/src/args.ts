export interface AncestralArgs {
  tree: string;
  outdir: string;
  input_fastas?: string[];
  aln?: string;
  vcf_reference?: string;
  alphabet?: string;
  model_name?: string;
  gtr_params?: string[];
  method_anc?: string;
  dense?: boolean;
  aa?: boolean;
  gap_fill?: string;
  keep_overhangs?: boolean;
  zero_based?: boolean;
  reconstruct_tip_states?: boolean;
  report_ambiguous?: boolean;
  gtr_iterations?: number;
  site_specific_gtr?: boolean;
  seed?: number;
}

export interface ClockArgs {
  dates: string;
  outdir: string;
  aln?: string[];
  tree?: string;
  vcf_reference?: string;
  name_column?: string;
  date_column?: string;
  sequence_length?: number;
  gtr?: string;
  gtr_params?: string[];
  branch_length_mode?: string;
  method_anc?: string;
  clock_filter?: number;
  reroot?: string;
  keep_root?: boolean;
  prune_short?: boolean;
  tip_slack?: number;
  covariation?: boolean;
  allow_negative_rate?: boolean;
  seed?: number;
}

export interface TimetreeArgs {
  outdir: string;
  tree?: string;
  dates?: string;
  input_fastas?: string[];
  aln?: string;
  vcf_reference?: string;
  clock_rate?: number;
  max_iter?: number;
  seed?: number;
}

export interface MugrationArgs {
  attribute: string;
  states: string;
  outdir: string;
  tree?: string;
  weights?: string;
  name_column?: string;
  confidence?: string;
  pc?: number;
  missing_data?: string;
  missing_weights_threshold?: number;
  iterations?: number;
  sampling_bias_correction?: number;
}

export interface OptimizeArgs {
  tree: string;
  outdir: string;
  input_fastas?: string[];
  alphabet?: string;
  model_name?: string;
  dense?: boolean;
  max_iter?: number;
  dp?: number;
  damping?: number;
  branch_length_initial_guess?: string;
  opt_method?: string;
  no_indels?: boolean;
}

export interface PruneArgs {
  tree: string;
  outdir: string;
  input_fastas?: string[];
  alphabet?: string;
  prune_short?: number;
  prune_empty?: boolean;
  merge_shared_mutations?: boolean;
  prune_nodes_list?: string;
  prune_nodes_list_delimiter?: string;
  prune_nodes_list_file?: string;
  prune_nodes_list_file_delimiter?: string;
}
