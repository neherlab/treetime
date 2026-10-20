use crate::make_error;
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use treetime_utils::fmt::string::counted;

const MIN_STATES: usize = 2;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StateAlphabet {
  states: IndexSet<String>,
}

impl StateAlphabet {
  pub fn new<S: Into<String>>(names: impl IntoIterator<Item = S>) -> Result<Self, Report> {
    let states: IndexSet<String> = names.into_iter().map(Into::into).sorted().dedup().collect();
    if states.len() < MIN_STATES {
      let found = if states.is_empty() {
        String::new()
      } else {
        format!(": {}", states.iter().join(", "))
      };
      return make_error!(
        "A discrete trait needs at least {MIN_STATES} different states, found {}{found}. \
         Check that the trait column holds at least {MIN_STATES} different known values.",
        counted(states.len(), "state")
      );
    }
    Ok(Self { states })
  }

  #[expect(clippy::len_without_is_empty, reason = "the constructor rejects fewer than 2 states")]
  pub fn len(&self) -> usize {
    self.states.len()
  }

  pub fn index_of(&self, name: &str) -> Option<usize> {
    self.states.get_index_of(name)
  }

  pub fn name(&self, index: usize) -> Option<&str> {
    self.states.get_index(index).map(String::as_str)
  }

  pub fn iter(&self) -> impl Iterator<Item = &str> {
    self.states.iter().map(String::as_str)
  }
}
