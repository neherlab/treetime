#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{
    AlignmentModel, GtrModelName, GtrOutput, JC69Params, ModelEntry, ModelOutput, SequenceModel, TraitModel,
    gtr_markdown, jc69,
  };
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_value_read_str, json_write_str};
  use treetime_utils::vec_of_owned;

  #[test]
  fn test_gtr_output_writes_the_model_fields() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let output = GtrOutput::builder().gtr(&gtr).model_name(GtrModelName::JC69).build();
    let parsed = json_value_read_str(&json_write_str(&output, JsonPretty(false))?)?;

    let fields = parsed
      .as_object()
      .expect("GTR JSON is an object")
      .keys()
      .cloned()
      .collect::<Vec<_>>();
    assert_eq!(vec!["model_type", "model_name", "mu", "pi", "W", "n_states"], fields);
    assert_eq!(parsed["n_states"], 4);
    assert_eq!(parsed["model_name"], "jc69");
    assert_eq!(parsed["model_type"], "named");
    Ok(())
  }

  #[test]
  fn test_gtr_markdown_writes_the_frequencies_and_rates_as_tables_with_the_states() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;

    let actual = gtr_markdown(&gtr, GtrModelName::JC69, &vec_of_owned!["A", "C", "G", "T"]);

    let expected = indoc! {"
      GTR model initialized: JC69, mu = 1.0000

      | pi | A | C | G | T |
      |---|---:|---:|---:|---:|
      | | 0.2500 | 0.2500 | 0.2500 | 0.2500 |

      | W | A | C | G | T |
      |---|---:|---:|---:|---:|
      | A | 0.0000 | 1.3333 | 1.3333 | 1.3333 |
      | C | 1.3333 | 0.0000 | 1.3333 | 1.3333 |
      | G | 1.3333 | 1.3333 | 0.0000 | 1.3333 |
      | T | 1.3333 | 1.3333 | 1.3333 | 0.0000 |"};
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_model_output_round_trips_sequence_and_trait_entries() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let model = GtrOutput::builder().gtr(&gtr).model_name(GtrModelName::Infer).build();
    let original = ModelOutput {
      alignments: vec![AlignmentModel {
        alignment: 1,
        path: PathBuf::from("na/aligned.fasta"),
        sites: 1410,
        indel_rate: Some(0.5),
      }],
      models: vec![
        ModelEntry::Sequence(SequenceModel {
          partition: 2,
          name: Some("cp3".to_owned()),
          alignment: 1,
          sites: 470,
          clock_rate: None,
          states: vec_of_owned!["A", "C", "G", "T"],
          model: model.clone(),
        }),
        ModelEntry::Trait(TraitModel {
          attribute: "region".to_owned(),
          states: vec_of_owned!["asia", "europe"],
          model,
        }),
      ],
    };

    let json = json_write_str(&original, JsonPretty(true))?;
    let restored: ModelOutput = json_read_str(&json)?;

    assert_eq!(json, json_write_str(&restored, JsonPretty(true))?);
    Ok(())
  }
}
