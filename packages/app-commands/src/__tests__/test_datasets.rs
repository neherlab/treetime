#[cfg(test)]
mod tests {
  use crate::check_inputs::InputKind;
  use crate::command::AppCommand;
  use crate::datasets::{DatasetInput, dataset_catalog};
  use pretty_assertions::assert_eq;
  use std::path::Path;
  use treetime_utils::o;

  #[test]
  fn test_datasets_zika_86_inputs_by_kind() {
    let catalog = dataset_catalog(Path::new("../../data")).unwrap();
    let zika = catalog
      .datasets
      .iter()
      .find(|dataset| dataset.name == "zika/86")
      .unwrap();
    assert_eq!(
      vec![
        DatasetInput {
          kind: InputKind::Tree,
          file: o!("zika/86/tree.nwk"),
          path: o!("../../data/zika/86/tree.nwk"),
        },
        DatasetInput {
          kind: InputKind::Alignment,
          file: o!("zika/86/aln.fasta.xz"),
          path: o!("../../data/zika/86/aln.fasta.xz"),
        },
        DatasetInput {
          kind: InputKind::Metadata,
          file: o!("zika/86/metadata.tsv"),
          path: o!("../../data/zika/86/metadata.tsv"),
        },
      ],
      zika.inputs
    );
  }

  #[test]
  fn test_datasets_examples_name_their_command() {
    let catalog = dataset_catalog(Path::new("../../data")).unwrap();
    let example = catalog
      .examples
      .iter()
      .find(|example| example.path == "ebola/20/ancestral-parsimony.yaml")
      .unwrap();
    assert_eq!(AppCommand::Ancestral, example.command);
  }
}
