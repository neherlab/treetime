#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::RecordingLog;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::sequence_inputs::{SourceSettings, read_sequence_sources};
  use eyre::Report;
  use helpers::{SourceFacts, facts, read_sources, read_sources_logged};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::path::PathBuf;
  use tempfile::tempdir;
  use treetime::ancestral::attach::MissingLeafLimit;
  use treetime::progress::LogLevel;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_sequence_sources_read_one_source_per_file_in_order() -> Result<(), Report> {
    let dir = tempdir()?;
    let first = dir.path().join("first.fasta");
    let second = dir.path().join("second.fasta");
    fs::write(&first, ">A one\nACGT\n>B\nACGA\n>C\nACGG\n>X\nAAAA\n")?;
    fs::write(&second, ">C\nTTN\n>A two\nTAN\n")?;

    let sources = read_sources(&[first.clone(), second.clone()], MissingLeafLimit::Enforced)?;

    let expected = vec![
      SourceFacts {
        path: first,
        length: 4,
        n_records: 4,
        mask: vec![false, false, false, false],
        seqs: btreemap! { o!("A") => o!("ACGT"), o!("B") => o!("ACGA"), o!("C") => o!("ACGG") },
        descs: btreemap! { o!("A") => Some(o!("one")), o!("B") => None, o!("C") => None },
      },
      SourceFacts {
        path: second,
        length: 3,
        n_records: 2,
        mask: vec![false, false, true],
        seqs: btreemap! { o!("A") => o!("TAN"), o!("B") => o!("NNN"), o!("C") => o!("TTN") },
        descs: btreemap! { o!("A") => Some(o!("two")), o!("C") => None },
      },
    ];
    assert_eq!(expected, facts(&sources));
    Ok(())
  }

  #[test]
  fn test_sequence_sources_apply_the_one_third_rule_per_file() -> Result<(), Report> {
    let dir = tempdir()?;
    let first = dir.path().join("first.fasta");
    let second = dir.path().join("second.fasta");
    fs::write(&first, ">A\nACGT\n>B\nACGA\n>C\nACGG\n")?;
    fs::write(&second, ">A\nTT\n")?;

    assert_error!(
      read_sources(
        &[first, second.clone()],
        MissingLeafLimit::Overridable { ignored: false }
      ),
      format!(
        "In alignment '{}': More than one third of terminal nodes (2 of 3) cannot be assigned a sequence. \
         Are you sure the alignment belongs to the tree? \
         Pass --allow-missing-sequences to proceed, treating missing tips as fully ambiguous.",
        second.display()
      )
    );
    Ok(())
  }

  #[test]
  fn test_sequence_sources_name_the_file_in_missing_leaf_warnings() -> Result<(), Report> {
    let dir = tempdir()?;
    let first = dir.path().join("first.fasta");
    let second = dir.path().join("second.fasta");
    fs::write(&first, ">A\nACGT\n>B\nACGA\n>C\nACGG\n")?;
    fs::write(&second, ">A\nTT\n>C\nTA\n")?;
    let log = RecordingLog::default();

    read_sources_logged(&[first, second.clone()], MissingLeafLimit::Enforced, &log)?;

    assert_eq!(
      vec![
        format!(
          "No sequence found for leaf 'B' in alignment '{}'; treating it as fully ambiguous (missing data).",
          second.display()
        ),
        format!(
          "1 of 3 tips have no matching sequence in alignment '{}' and are treated as fully ambiguous.",
          second.display()
        ),
      ],
      log.messages(LogLevel::Warn)
    );
    Ok(())
  }

  #[test]
  fn test_sequence_sources_reject_standard_input_given_twice() {
    assert_error!(
      read_sources(&[PathBuf::from("-"), PathBuf::from("-")], MissingLeafLimit::Enforced),
      "Standard input ('-') can be read only once, but --alignment lists it 2 times"
    );
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::sequence_inputs::AlignmentSource;
    use std::collections::BTreeMap;
    use treetime::alphabet::alphabet::Alphabet;
    use treetime::progress::{LogSink, NoopProgress};
    use treetime_io::nwk::nwk_read;

    #[derive(Debug, PartialEq)]
    pub(super) struct SourceFacts {
      pub(super) path: PathBuf,
      pub(super) length: usize,
      pub(super) n_records: usize,
      pub(super) mask: Vec<bool>,
      pub(super) seqs: BTreeMap<String, String>,
      pub(super) descs: BTreeMap<String, Option<String>>,
    }

    pub(super) fn read_sources(paths: &[PathBuf], limit: MissingLeafLimit) -> Result<Vec<AlignmentSource>, Report> {
      read_sources_logged(paths, limit, &NoopProgress)
    }

    pub(super) fn read_sources_logged(
      paths: &[PathBuf],
      limit: MissingLeafLimit,
      log: &dyn LogSink,
    ) -> Result<Vec<AlignmentSource>, Report> {
      let tree = nwk_read(b"(A:0.1,B:0.1,C:0.1)root;".as_slice())?;
      let names = tree.names.clone();
      let alphabet = Alphabet::default();
      let settings = SourceSettings {
        alphabet: &alphabet,
        gap_fill: None,
        missing_leaves: limit,
      };
      let alignment = AlignmentArgs {
        alignment: paths.to_vec(),
      };
      read_sequence_sources(&alignment, &settings, &tree.graph, &names, log)
    }

    pub(super) fn facts(sources: &[AlignmentSource]) -> Vec<SourceFacts> {
      sources
        .iter()
        .map(|aligned| {
          let source = &aligned.source;
          let name = |key| source.nodes[key].name.clone();
          SourceFacts {
            path: source.path.clone(),
            length: source.length,
            n_records: aligned.n_records,
            mask: source.mask.clone(),
            seqs: source
              .nodes
              .values()
              .filter_map(|node| Some((node.name.clone(), node.seq.as_ref()?.to_string())))
              .collect(),
            descs: aligned
              .descs
              .iter()
              .map(|(key, desc)| (name(key), desc.clone()))
              .collect(),
          }
        })
        .collect()
    }
  }
}
