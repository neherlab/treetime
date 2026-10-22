#[cfg(test)]
mod tests {
  use crate::command::CommandArgs;
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
  use crate::commands::shared::tree_input::TreeArgs;
  use std::path::PathBuf;
  use treetime_utils::assert_error;

  #[test]
  fn test_alignment_input_is_required() {
    let result = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
      tree: TreeArgs {
        tree: Some(PathBuf::from("tree.nwk")),
      },
      reconstructed_fasta_output: ReconstructedFastaOutputArgs {
        output_sequences: Some(PathBuf::from("out.fasta")),
      },
      ..TreetimeAncestralArgsRaw::default()
    });

    assert_error!(
      result,
      "--alignment is required: pass one or more FASTA files, or '-' to read standard input"
    );
  }

  #[test]
  fn test_alignment_input_reads_standard_input_once() {
    let raw = TreetimeAncestralArgsRaw {
      tree: TreeArgs {
        tree: Some(PathBuf::from("tree.nwk")),
      },
      alignment: AlignmentArgs {
        alignment: vec![PathBuf::from("-"), PathBuf::from("-")],
      },
      reconstructed_fasta_output: ReconstructedFastaOutputArgs {
        output_sequences: Some(PathBuf::from("out.{alignment}.fasta")),
      },
      ..TreetimeAncestralArgsRaw::default()
    };

    assert_error!(
      CommandArgs::try_from(raw),
      "Standard input can be read only once, but --alignment lists it 2 times"
    );
  }
}
