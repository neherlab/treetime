use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use treetime_primitives::AsciiChar;
use treetime_utils::make_error;

/// Gene key for nucleotide substitutions.
///
/// Matches the key augur uses in `branch_attrs.mutations` and in the
/// `root_sequence` map for the whole-genome nucleotide track.
pub const NUC_GENE: &str = "nuc";

/// A single substitution on a branch: ancestral state -> derived state at one position.
///
/// `gene` is [`NUC_GENE`] for nucleotide substitutions or a CDS/gene name for
/// amino-acid substitutions. `position` is 1-based, matching both the Auspice
/// mutation-string convention (`A123T`) and the UShER protobuf position field.
/// `parent` and `child` are the ancestral and derived characters (`A`/`C`/`G`/`T`
/// for nucleotides, single-letter amino acids otherwise).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct TreeIrSub {
  pub gene: String,
  pub position: usize,
  pub parent: AsciiChar,
  pub child: AsciiChar,
}

impl TreeIrSub {
  pub fn is_nuc(&self) -> bool {
    self.gene == NUC_GENE
  }

  /// Render in Auspice mutation-string form, e.g. `A123T`.
  pub fn to_auspice_string(&self) -> String {
    format!("{}{}{}", self.parent, self.position, self.child)
  }

  /// Parse an Auspice mutation string (`A123T`) for a given gene track.
  ///
  /// The leading character is the ancestral state, the trailing character is the
  /// derived state, and the digits in between are the 1-based position.
  pub fn from_auspice_string(gene: &str, s: &str) -> Result<Self, Report> {
    let bytes = s.as_bytes();
    if bytes.len() < 3 {
      return make_error!("Invalid mutation string '{s}': expected form like 'A123T'");
    }
    let parent = AsciiChar::try_new(bytes[0])?;
    let child = AsciiChar::try_new(bytes[bytes.len() - 1])?;
    let position_bytes = &bytes[1..bytes.len() - 1];
    let position: usize = std::str::from_utf8(position_bytes)
      .wrap_err_with(|| format!("Invalid mutation string '{s}': position is not valid UTF-8"))?
      .parse()
      .wrap_err_with(|| format!("Invalid mutation string '{s}': position is not a number"))?;
    Ok(Self {
      gene: gene.to_owned(),
      position,
      parent,
      child,
    })
  }
}

/// Whether an indel is an insertion or a deletion relative to the parent.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum IndelKind {
  Insertion,
  Deletion,
}

/// An insertion or deletion event on a branch.
///
/// `start` is the 1-based start position; `seq` is the inserted (for
/// [`IndelKind::Insertion`]) or deleted (for [`IndelKind::Deletion`]) sequence.
/// Indels are not representable in every output format (UShER MAT carries only
/// nucleotide substitutions); adapters that cannot represent them drop indels
/// with a warning.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct TreeIrIndel {
  pub gene: String,
  pub kind: IndelKind,
  pub start: usize,
  pub seq: Vec<AsciiChar>,
}
