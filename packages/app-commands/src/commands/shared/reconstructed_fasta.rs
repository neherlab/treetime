#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaWriter;
use treetime_primitives::Seq;
use treetime_schema::{schema_defaults, skip_serializing_optionals};
use treetime_utils::make_internal_error;

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ReconstructedFastaOutputArgs {
  /// Path to output reconstructed nucleotide FASTA.
  ///
  /// Holds the reconstructed sequences of the internal nodes, and of the tips with `--include-leaves` or
  /// `--reconstruct-tip-states`; `--impute-missing-data` resolves ambiguous tip states. Takes precedence over paths
  /// configured with `--output-all` and `--output-selection`.
  ///
  /// With several `--alignment` files, one file is written per alignment, and the file name must contain
  /// `{partition}`, which is replaced by the number of the alignment (1, 2, ...).
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Output"))]
  #[schemars(extend("x-path" = "output"))]
  pub output_reconstructed_nuc_fasta: Option<PathBuf>,
}

pub(crate) struct NucFastaWriters {
  partitions: BTreeMap<usize, NucFasta>,
}

impl NucFastaWriters {
  pub(crate) fn create(
    paths: BTreeMap<usize, PathBuf>,
    mut descs: BTreeMap<usize, BTreeMap<GraphNodeKey, Option<String>>>,
  ) -> Result<Self, Report> {
    let partitions = paths
      .into_iter()
      .map(|(partition, path)| {
        let fasta = NucFasta {
          writer: FastaWriter::create(&path)?,
          descs: descs.remove(&partition).unwrap_or_default(),
          path,
        };
        Ok((partition, fasta))
      })
      .collect::<Result<_, Report>>()?;
    Ok(Self { partitions })
  }

  pub(crate) fn write(&mut self, partition: usize, key: GraphNodeKey, name: &str, seq: &Seq) -> Result<(), Report> {
    let Some(fasta) = self.partitions.get_mut(&partition) else {
      return make_internal_error!("No reconstructed nucleotide FASTA was planned for partition {partition}");
    };
    let desc = fasta.descs.get(&key).and_then(Option::as_deref);
    fasta.writer.write(name, desc, seq)
  }

  pub(crate) fn finish(self) -> Result<Vec<PathBuf>, Report> {
    self
      .partitions
      .into_values()
      .map(|fasta| {
        fasta.writer.finish()?;
        Ok(fasta.path)
      })
      .collect()
  }
}

struct NucFasta {
  path: PathBuf,
  writer: FastaWriter,
  descs: BTreeMap<GraphNodeKey, Option<String>>,
}
