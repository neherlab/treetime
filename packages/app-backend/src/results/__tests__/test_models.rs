#[cfg(test)]
mod tests {
  use crate::results::models::{RunModel, run_models};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime::gtr::get_gtr::ModelOutput;
  use treetime_utils::io::json::json_read_str;

  #[test]
  fn test_run_models_give_the_frequencies_and_the_rates_from_each_state_to_each_other() -> Result<(), Report> {
    let output: ModelOutput = json_read_str(indoc! {r#"{
      "alignments": [{ "alignment": 1, "path": "aln.fasta", "sites": 100 }],
      "models": [{
        "kind": "sequence",
        "partition": 1,
        "alignment": 1,
        "sites": 100,
        "states": ["A", "C", "G"],
        "model_type": "inferred",
        "model_name": "infer",
        "mu": 2.0,
        "pi": [0.5, 0.25, 0.25],
        "W": [[0.0, 1.0, 4.0], [1.0, 0.0, 2.0], [4.0, 2.0, 0.0]],
        "n_states": 3
      }]
    }"#})?;

    let models = run_models(Some(&output));

    let RunModel::Sequence(model) = &models.models[0] else {
      panic!("expected a sequence model");
    };
    assert_eq!(
      (
        vec![0.5, 0.25, 0.25],
        vec![vec![0.0, 0.5, 2.0], vec![1.0, 0.0, 1.0], vec![4.0, 1.0, 0.0]],
      ),
      (model.pi.clone(), model.rates.clone())
    );
    Ok(())
  }
}
