#[cfg(test)]
mod tests {
  use crate::commands::shared::gtr_output::{trait_model, write_gtr_output};
  use app_output::output_plan::Requested;
  use eyre::Report;
  use helpers::{gtr_output, gtr_with_mu, jc69_gtr, jc69_w, models, partition, read_document, split_partitions};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use tempfile::TempDir;
  use treetime::progress::NoopProgress;
  use treetime_utils::{assert_error, vec_of_owned};

  #[test]
  fn test_gtr_output_one_partition_writes_one_sequence_entry() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let resolved = gtr_output(&dir, Requested::Named);
    let results = vec![partition(1, "aln.fasta", Some(jc69_gtr()?), None)];

    write_gtr_output(&resolved, models(results, None)?, "no model", &NoopProgress)?;

    let expected = json!({
      "models": [{
        "kind": "sequence",
        "partition": 1,
        "alignment": "aln.fasta",
        "sites": 12,
        "model_type": "named",
        "model_name": "jc69",
        "mu": 1.0,
        "pi": [0.25, 0.25, 0.25, 0.25],
        "W": jc69_w(),
        "n_states": 4,
      }],
    });
    assert_eq!(expected, read_document(&resolved)?);
    Ok(())
  }

  #[test]
  fn test_gtr_output_two_partitions_write_entries_in_partition_order_with_file_clock_rates() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let resolved = gtr_output(&dir, Requested::Named);
    let results = vec![
      partition(1, "ha/aligned.fasta", Some(gtr_with_mu(0.5)?), Some(0.25)),
      partition(2, "na/aligned.fasta", Some(gtr_with_mu(2.0)?), Some(0.125)),
    ];

    write_gtr_output(&resolved, models(results, Some(0.004))?, "no model", &NoopProgress)?;

    let document = read_document(&resolved)?;
    let summary = document["models"]
      .as_array()
      .expect("models is an array")
      .iter()
      .map(|entry| {
        json!([
          entry["partition"],
          entry["alignment"],
          entry["mu"],
          entry["indel_rate"],
          entry["clock_rate"]
        ])
      })
      .collect::<Vec<_>>();
    assert_eq!(
      vec![
        json!([1, "ha/aligned.fasta", 0.5, 0.25, 0.002]),
        json!([2, "na/aligned.fasta", 2.0, 0.125, 0.008]),
      ],
      summary
    );
    Ok(())
  }

  #[test]
  fn test_gtr_output_partitions_of_one_alignment_take_its_path_and_indel_rate() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let resolved = gtr_output(&dir, Requested::Named);

    write_gtr_output(&resolved, split_partitions()?, "no model", &NoopProgress)?;

    let summary = read_document(&resolved)?["models"]
      .as_array()
      .expect("models is an array")
      .iter()
      .map(|entry| {
        json!([
          entry["partition"],
          entry["alignment"],
          entry["sites"],
          entry["indel_rate"]
        ])
      })
      .collect::<Vec<_>>();
    assert_eq!(
      vec![json!([1, "genome.fasta", 7, 0.5]), json!([2, "genome.fasta", 5, 0.5])],
      summary
    );
    Ok(())
  }

  #[test]
  fn test_gtr_output_fitch_partition_has_no_entry() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let resolved = gtr_output(&dir, Requested::Named);
    let results = vec![
      partition(1, "a.fasta", None, None),
      partition(2, "b.fasta", Some(jc69_gtr()?), None),
    ];

    write_gtr_output(&resolved, models(results, None)?, "no model", &NoopProgress)?;

    let partitions = read_document(&resolved)?["models"]
      .as_array()
      .expect("models is an array")
      .iter()
      .map(|entry| entry["partition"].clone())
      .collect::<Vec<_>>();
    assert_eq!(vec![json!(2)], partitions);
    Ok(())
  }

  #[test]
  fn test_gtr_output_fitch_only_run_fails_when_the_file_is_named() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let resolved = gtr_output(&dir, Requested::Named);
    let results = vec![partition(1, "a.fasta", None, None)];

    assert_error!(
      write_gtr_output(
        &resolved,
        models(results, None)?,
        "no GTR model was fitted",
        &NoopProgress
      ),
      "--output-model was requested, but no GTR model was fitted"
    );
    Ok(())
  }

  #[test]
  fn test_gtr_output_fitch_only_run_skips_the_file_of_output_all() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let resolved = gtr_output(&dir, Requested::All);
    let results = vec![partition(1, "a.fasta", None, None)];

    write_gtr_output(&resolved, models(results, None)?, "no model", &NoopProgress)?;

    assert!(!dir.path().join("model.json").exists());
    Ok(())
  }

  #[test]
  fn test_gtr_output_mugration_entry_keeps_attribute_and_states() -> Result<(), Report> {
    let dir = TempDir::new()?;
    let resolved = gtr_output(&dir, Requested::Named);
    let gtr = jc69_gtr()?;

    write_gtr_output(
      &resolved,
      vec![trait_model(&gtr, "country", vec_of_owned!["germany", "usa"])],
      "no model",
      &NoopProgress,
    )?;

    let expected = json!({
      "models": [{
        "kind": "trait",
        "attribute": "country",
        "states": ["germany", "usa"],
        "model_type": "inferred",
        "model_name": "infer",
        "mu": 1.0,
        "pi": [0.25, 0.25, 0.25, 0.25],
        "W": jc69_w(),
        "n_states": 4,
      }],
    });
    assert_eq!(expected, read_document(&resolved)?);
    Ok(())
  }

  mod helpers {
    use crate::commands::shared::gtr_output::partition_models;
    use app_output::output_plan::{OutputSelection, PlannedFile, Requested, ResolvedOutputs};
    use eyre::Report;
    use maplit::btreemap;
    use serde_json::{Value, json};
    use std::num::NonZeroUsize;
    use std::path::PathBuf;
    use tempfile::TempDir;
    use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
    use treetime::gtr::get_gtr::{GtrModelName, JC69Params, ModelEntry, jc69};
    use treetime::gtr::gtr::{GTR, Mu};
    use treetime::partition::result::{AlignmentResult, PartitionResult};
    use treetime_utils::io::json::json_value_read_file;

    pub(super) fn gtr_output(dir: &TempDir, requested: Requested) -> ResolvedOutputs {
      ResolvedOutputs {
        alignment_count: NonZeroUsize::MIN,
        tree_outputs: btreemap! {},
        non_tree_outputs: btreemap! {
          OutputSelection::Model => PlannedFile {
            path: dir.path().join("model.json"),
            requested,
            condition: None,
          },
        },
        tree_outputs_by_alignment: btreemap! {},
        non_tree_outputs_by_alignment: btreemap! {},
      }
    }

    pub(super) fn read_document(resolved: &ResolvedOutputs) -> Result<Value, Report> {
      json_value_read_file(&resolved.non_tree_outputs[&OutputSelection::Model].path)
    }

    pub(super) fn jc69_gtr() -> Result<GTR, Report> {
      jc69(JC69Params::default())
    }

    pub(super) fn jc69_w() -> Value {
      let w = 4.0 / 3.0;
      json!([[0.0, w, w, w], [w, 0.0, w, w], [w, w, 0.0, w], [w, w, w, 0.0]])
    }

    pub(super) fn gtr_with_mu(mu: f64) -> Result<GTR, Report> {
      let base = jc69_gtr()?;
      GTR::builder()
        .n_states(4)
        .mu(Mu::Normalized(mu))
        .W(base.W().clone())
        .pi(base.pi().clone())
        .build()
    }

    pub(super) fn partition(
      number: usize,
      path: &str,
      gtr: Option<GTR>,
      indel_rate: Option<f64>,
    ) -> (PartitionResult, AlignmentResult) {
      let partition = PartitionResult {
        number,
        name: None,
        alignment: number,
        sites: 12,
        model_name: if gtr.is_some() {
          GtrModelName::JC69
        } else {
          GtrModelName::Infer
        },
        gtr,
        log_lh: None,
      };
      (partition, alignment(number, path, indel_rate))
    }

    pub(super) fn models(
      results: Vec<(PartitionResult, AlignmentResult)>,
      clock_rate: Option<f64>,
    ) -> Result<Vec<ModelEntry>, Report> {
      let (partitions, alignments): (Vec<_>, Vec<_>) = results.into_iter().unzip();
      partition_models(&partitions, &alignments, clock_rate)
    }

    pub(super) fn split_partitions() -> Result<Vec<ModelEntry>, Report> {
      let partition = |number, sites| -> Result<PartitionResult, Report> {
        Ok(PartitionResult {
          number,
          name: None,
          alignment: 1,
          sites,
          model_name: GtrModelName::JC69,
          gtr: Some(jc69_gtr()?),
          log_lh: None,
        })
      };
      partition_models(
        &[partition(1, 7)?, partition(2, 5)?],
        &[alignment(1, "genome.fasta", Some(0.5))],
        None,
      )
    }

    fn alignment(number: usize, path: &str, indel_rate: Option<f64>) -> AlignmentResult {
      AlignmentResult {
        number,
        path: PathBuf::from(path),
        alphabet: Alphabet::new(AlphabetName::Nuc).unwrap(),
        length: 12,
        mask: vec![true; 12],
        indel_rate,
        sequences: None,
      }
    }
  }
}
