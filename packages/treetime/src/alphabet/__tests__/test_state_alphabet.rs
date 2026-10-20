#[cfg(test)]
mod tests {
  use crate::alphabet::state_alphabet::StateAlphabet;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use treetime_utils::assert_error;

  #[test]
  fn test_state_alphabet_sorts_and_deduplicates() -> Result<(), Report> {
    let states = StateAlphabet::new(["USA", "Germany", "China", "USA"])?;

    assert_eq!(3, states.len());
    assert_eq!(vec!["China", "Germany", "USA"], states.iter().collect_vec());
    Ok(())
  }

  #[test]
  fn test_state_alphabet_index_of_and_name_agree() -> Result<(), Report> {
    let states = StateAlphabet::new(["USA", "Germany", "China"])?;

    let expected = vec![
      (Some(0), Some("China")),
      (Some(1), Some("Germany")),
      (Some(2), Some("USA")),
    ];
    let actual = ["China", "Germany", "USA"]
      .iter()
      .map(|name| {
        let index = states.index_of(name);
        (index, index.and_then(|index| states.name(index)))
      })
      .collect_vec();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_state_alphabet_unknown_name_and_index() -> Result<(), Report> {
    let states = StateAlphabet::new(["USA", "Germany"])?;

    assert_eq!(None, states.index_of("France"));
    assert_eq!(None, states.name(2));
    Ok(())
  }

  #[test]
  fn test_state_alphabet_rejects_zero_states() {
    assert_error!(
      StateAlphabet::new(Vec::<String>::new()),
      "A discrete trait needs at least 2 different states, found 0 states. \
       Check that the trait column holds at least 2 different known values."
    );
  }

  #[test]
  fn test_state_alphabet_rejects_one_state() {
    assert_error!(
      StateAlphabet::new(["usa", "usa"]),
      "A discrete trait needs at least 2 different states, found 1 state: usa. \
       Check that the trait column holds at least 2 different known values."
    );
  }
}
