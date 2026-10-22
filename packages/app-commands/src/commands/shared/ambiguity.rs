use crate::commands::shared::argument_error::argument_error;
use app_output::ambiguity_outputs::AmbiguityOutputs;
use app_output::output_plan::OutputSelection;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use std::collections::BTreeSet;
use strum::IntoEnumIterator;
use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

pub fn resolve_ambiguity_outputs<O>(
  ambiguity_threshold: Option<AmbiguityThreshold>,
  ambiguity_outputs: &[O],
) -> Result<AmbiguityOutputs, Report>
where
  O: Copy + Into<OutputSelection> + IntoEnumIterator,
{
  let selected: BTreeSet<OutputSelection> = ambiguity_outputs.iter().copied().map(Into::into).collect();
  match ambiguity_threshold {
    Some(threshold) => Ok(AmbiguityOutputs::new(threshold, &selected, O::iter().map(Into::into))),
    None if selected.is_empty() => Ok(AmbiguityOutputs::default()),
    None => Err(argument_error(
      &["ambiguity_threshold", "ambiguity_outputs"],
      "--ambiguity-outputs requires --ambiguity-threshold",
    )),
  }
}

macro_rules! ambiguity_args {
  ($args:ident, $output:ident, $details:literal { $($(#[$doc:meta])* $variant:ident),* $(,)? }) => {
    #[derive(Copy, Clone, Debug, Eq, PartialEq, JsonSchema, strum_macros::EnumIter, Serialize, Deserialize)]
    #[schemars(rename_all = "kebab-case")]
    #[deser(rename_all = "kebab-case")]
    #[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
    pub enum $output {
      $($(#[$doc])* $variant),*
    }

    impl From<$output> for OutputSelection {
      fn from(value: $output) -> Self {
        match value {
          $($output::$variant => Self::$variant),*
        }
      }
    }

    #[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
    #[deser(skip_serializing_optionals)]
    #[schemars(transform = skip_serializing_optionals)]
    #[schemars(default, deny_unknown_fields)]
    #[schemars(transform = schema_defaults::<Self>)]
    #[deser(default, deny_unknown_fields)]
    #[cfg_attr(feature = "clap", derive(clap::Args))]
    pub struct $args {
      /// Report every state whose posterior probability is at most this far below the probability of the most likely
      /// state.
      ///
      /// The value is a difference of probabilities from 0 to 1. With 0, the report names only the states that tie with
      /// the most likely state. RAxML-NG uses 0.5 divided by the number of states, for example 0.125 for nucleotides.
      /// The report is off when the value is not set.
      ///
      #[doc = $details]
      #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, help_heading = "Output"))]
      pub ambiguity_threshold: Option<AmbiguityThreshold>,

      /// List of the outputs that carry the ambiguity report.
      ///
      /// Without this list, every output of the run that can carry the report carries it. Requires
      /// `--ambiguity-threshold`.
      #[cfg_attr(feature = "clap", clap(long, hide_short_help = true, num_args = 1.., value_delimiter = ',', help_heading = "Output"))]
      pub ambiguity_outputs: Vec<$output>,
    }
  };
}

ambiguity_args!(AncestralAmbiguityArgs, AncestralAmbiguityOutput,
  "A sequence output writes a position with more than one reported state as the ambiguity code of these states: an \
   IUPAC code for nucleotides, `B`, `Z`, or `J` for amino acids, otherwise `N` or `X`. Internal nodes have a report \
   at every position, leaves only at the positions that `--impute-missing-data` imputes. Nodes whose sequence \
   `--sample-from-profile` draws have no report. The setting applies only when a marginal reconstruction runs \
   (`--method-anc marginal`, or the amino-acid reconstruction of `--translations`), because Fitch parsimony \
   computes no posterior." {
  /// Reconstructed nucleotide FASTA
  Sequences,
  /// Reconstructed amino-acid FASTA of each CDS
  Translations,
  /// Node sequences in the augur node data
  NodeData,
});

ambiguity_args!(TimetreeAmbiguityArgs, TimetreeAmbiguityOutput,
  "The reconstructed FASTA writes a position with more than one reported state as the IUPAC code of these states, \
   or `N` for all four nucleotides. Internal nodes have a report at every position, leaves only at the positions \
   that `--impute-missing-data` imputes." {
  /// Reconstructed nucleotide FASTA
  Sequences,
});

ambiguity_args!(MugrationAmbiguityArgs, MugrationAmbiguityOutput,
  "A trait output names the reported states of a node with more than one reported state, joined by `|` in the \
   order of the state names, for example `American Samoa|Samoa`. A node with one reported state has an empty cell in \
   the traits table and no field in the JSON outputs." {
  /// Extra column `<attribute>_ambiguous` in the traits table
  Traits,
  /// Extra node field `<attribute>_ambiguous` in the augur node data
  NodeData,
  /// Extra node attribute `<attribute>_ambiguous` in the Auspice JSON
  Auspice,
});
