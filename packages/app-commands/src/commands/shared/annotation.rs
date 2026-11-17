use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use std::path::PathBuf;
use treetime_schema::schema_defaults;

#[cfg(feature = "clap")]
use clap::ValueHint;

/// Genome annotations of the `--alignment` files
#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct AnnotationArgs {
  /// Genome annotation (GFF3) with the coordinates of the CDSs of an `--alignment` file
  ///
  /// With `--partition-by codon-positions` or `genes`, give one file per `--alignment` file, matched in list order:
  /// the scheme builds the partitions of each alignment from the CDSs of its annotation. In `ancestral`, the file also
  /// supplies the annotations of the augur node data and, without `--cds-selection`, the CDSs that `--translations`
  /// reconstructs; `--translations` takes one `--alignment` file and one annotation file. TreeTime reads the coding
  /// regions from the CDS rows, which need strand `+` or `-` and phase 0, 1 or 2. A path that ends in `.gz`, `.bz2`,
  /// `.xz`, or `.zst` reads a compressed file.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_hint = ValueHint::FilePath, hide_short_help = true, help_heading = "Input data")
  )]
  #[schemars(extend("x-path" = "input"))]
  pub annotation: Vec<PathBuf>,
}
