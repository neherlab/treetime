use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::gtr::get_gtr::{AlphabetRequirement, GtrModel};
use crate::partition::sites::SiteRange;
use derive_more::{Display, Error};
use eyre::Report;
use itertools::{Itertools, izip};
use treetime_utils::fmt::string::counted;
use treetime_utils::o;

pub fn resolve_partitions(
  entries: &[PartitionEntry],
  lengths: &[usize],
  default_model: &GtrModel,
  alphabet: Option<AlphabetName>,
) -> Result<ResolvedPartitions, Report> {
  let ResolvedLayouts {
    layouts,
    duplicate_names,
  } = resolve_partition_layouts(entries, lengths)?;
  if layouts.is_empty() {
    check_model_alphabet(default_model, model_alphabet(default_model, alphabet))?;
  }
  let alphabet = alphabet.unwrap_or_default();
  let n_states = Alphabet::new(alphabet)?.n_canonical();
  let plans = izip!(0.., layouts)
    .map(|(index, layout)| {
      let model = entries.get(index).map_or(default_model, |entry| &entry.model);
      if let Some(problem) = model_alphabet_problem(model, alphabet, n_states) {
        return model_error(problem.message(Some(&layout.label()), model, alphabet));
      }
      Ok(PartitionPlan {
        layout,
        model: model.clone(),
      })
    })
    .collect::<Result<_, Report>>()?;
  Ok(ResolvedPartitions { plans, duplicate_names })
}

pub fn model_alphabet(model: &GtrModel, alphabet: Option<AlphabetName>) -> AlphabetName {
  alphabet.or_else(|| model.alphabet()).unwrap_or_default()
}

pub fn check_model_alphabet(model: &GtrModel, alphabet: AlphabetName) -> Result<(), Report> {
  let n_states = Alphabet::new(alphabet)?.n_canonical();
  match model_alphabet_problem(model, alphabet, n_states) {
    Some(problem) => model_error(problem.message(None, model, alphabet)),
    None => Ok(()),
  }
}

pub fn resolve_partition_layouts(entries: &[PartitionEntry], lengths: &[usize]) -> Result<ResolvedLayouts, Report> {
  if entries.is_empty() {
    return whole_alignment_layouts(lengths);
  }
  if lengths.is_empty() {
    return layout_error(o!(
      "--partition needs --alignment: each entry takes sites of an --alignment file, but the run has no --alignment \
       file"
    ));
  }
  let alignments: Vec<usize> = izip!(entries, 1..)
    .map(|(entry, number)| entry_alignment(entry, number, lengths.len()))
    .try_collect()?;
  if let Some(alignment) = (1..=lengths.len()).find(|alignment| !alignments.contains(alignment)) {
    return layout_error(format!(
      "alignment {alignment} has no partition entry: with --partition, each --alignment file needs at least one \
       entry"
    ));
  }
  let mut owners: Vec<Vec<Option<usize>>> = lengths.iter().map(|&length| vec![None; length]).collect();
  for (index, (entry, &alignment)) in izip!(entries, &alignments).enumerate() {
    if let Some(sites) = &entry.sites {
      take_sites(entries, index, sites, alignment, &mut owners[alignment - 1])?;
    }
  }
  for (alignment, owners) in izip!(1.., &mut owners) {
    take_remaining_sites(entries, &alignments, alignment, owners)?;
  }
  let mut columns = vec![vec![]; entries.len()];
  for owners in &owners {
    for (column, owner) in owners.iter().enumerate() {
      if let Some(owner) = owner {
        columns[*owner].push(column);
      }
    }
  }
  let layouts = izip!(entries, alignments, columns, 1..)
    .map(|(entry, alignment, columns, number)| PartitionLayout {
      number,
      name: entry.name.clone(),
      alignment,
      columns,
    })
    .collect();
  Ok(ResolvedLayouts {
    layouts,
    duplicate_names: duplicate_names(entries),
  })
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedPartitions {
  pub plans: Vec<PartitionPlan>,
  pub duplicate_names: Vec<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PartitionPlan {
  pub layout: PartitionLayout,
  pub model: GtrModel,
}

fn whole_alignment_layouts(lengths: &[usize]) -> Result<ResolvedLayouts, Report> {
  let layouts = izip!(lengths, 1..)
    .map(|(&length, number)| {
      if length == 0 {
        return layout_error(format!("alignment {number} has no sites"));
      }
      Ok(PartitionLayout::whole(number, length))
    })
    .try_collect()?;
  Ok(ResolvedLayouts {
    layouts,
    duplicate_names: vec![],
  })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolvedLayouts {
  pub layouts: Vec<PartitionLayout>,
  pub duplicate_names: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PartitionLayout {
  pub number: usize,
  pub name: Option<String>,
  pub alignment: usize,
  pub columns: Vec<usize>,
}

impl PartitionLayout {
  pub fn whole(number: usize, length: usize) -> Self {
    Self {
      number,
      name: None,
      alignment: number,
      columns: (0..length).collect(),
    }
  }

  pub fn covers(&self, length: usize) -> bool {
    self.columns.len() == length
  }

  pub fn label(&self) -> String {
    partition_label(self.number, self.name.as_deref())
  }
}

fn entry_alignment(entry: &PartitionEntry, number: usize, n_alignments: usize) -> Result<usize, Report> {
  let label = || partition_label(number, entry.name.as_deref());
  match entry.alignment {
    Some(alignment) if (1..=n_alignments).contains(&alignment) => Ok(alignment),
    Some(alignment) => layout_error(format!(
      "{} names alignment {alignment}, but --alignment lists {}",
      label(),
      counted(n_alignments, "file")
    )),
    None if n_alignments == 1 => Ok(1),
    None => layout_error(format!(
      "{} needs `alignment`, the number of its --alignment file, because --alignment lists {}",
      label(),
      counted(n_alignments, "file")
    )),
  }
}

fn take_sites(
  entries: &[PartitionEntry],
  index: usize,
  sites: &[SiteRange],
  alignment: usize,
  owners: &mut [Option<usize>],
) -> Result<(), Report> {
  let label = |index: usize| partition_label(index + 1, entries[index].name.as_deref());
  if sites.is_empty() {
    return layout_error(format!(
      "{} has an empty `sites` list: list its sites, or leave out `sites` to take the remaining sites of its \
       alignment",
      label(index)
    ));
  }
  for range in sites {
    if range.last() > owners.len() {
      return layout_error(format!(
        "{}: the site range '{range}' ends after the last site {} of alignment {alignment}",
        label(index),
        owners.len()
      ));
    }
    for column in range.columns() {
      match owners[column] {
        Some(owner) if owner == index => {
          return layout_error(format!("{} lists site {} twice", label(index), column + 1));
        },
        Some(owner) => {
          return layout_error(format!(
            "{} takes site {}, which {} already takes",
            label(index),
            column + 1,
            label(owner)
          ));
        },
        None => owners[column] = Some(index),
      }
    }
  }
  Ok(())
}

fn take_remaining_sites(
  entries: &[PartitionEntry],
  alignments: &[usize],
  alignment: usize,
  owners: &mut [Option<usize>],
) -> Result<(), Report> {
  let label = |index: usize| partition_label(index + 1, entries[index].name.as_deref());
  let remainders = izip!(entries, alignments)
    .positions(|(entry, &entry_alignment)| entry_alignment == alignment && entry.sites.is_none())
    .collect_vec();
  let uncovered = owners.iter().positions(Option::is_none).collect_vec();
  match (remainders.as_slice(), uncovered.first()) {
    ([], None) => Ok(()),
    ([], Some(first)) => layout_error(format!(
      "no partition takes {} of alignment {alignment}, the first is site {}: list them in an entry, or add an entry \
       without `sites`, which takes the remaining sites",
      counted(uncovered.len(), "site"),
      first + 1
    )),
    ([index], None) => layout_error(format!(
      "{} leaves out `sites`, but the other entries take every site of alignment {alignment}",
      label(*index)
    )),
    ([index], Some(_)) => {
      for column in uncovered {
        owners[column] = Some(*index);
      }
      Ok(())
    },
    ([first, second, ..], _) => layout_error(format!(
      "{} and {} both leave out `sites` in alignment {alignment}, but only one entry of an alignment can take its \
       remaining sites",
      label(*first),
      label(*second)
    )),
  }
}

pub fn partition_label(number: usize, name: Option<&str>) -> String {
  match name {
    Some(name) => format!("partition {number} '{name}'"),
    None => format!("partition {number}"),
  }
}

fn duplicate_names(entries: &[PartitionEntry]) -> Vec<String> {
  entries
    .iter()
    .filter_map(|entry| entry.name.clone())
    .duplicates()
    .sorted()
    .collect()
}

#[derive(Clone, Debug, PartialEq)]
pub struct PartitionEntry {
  pub name: Option<String>,
  pub alignment: Option<usize>,
  pub sites: Option<Vec<SiteRange>>,
  pub model: GtrModel,
}

fn model_alphabet_problem(model: &GtrModel, alphabet: AlphabetName, n_states: usize) -> Option<ModelAlphabetProblem> {
  if let Some(requirement) = model
    .name()
    .alphabet_requirement()
    .filter(|requirement| requirement.alphabet != alphabet)
  {
    return Some(ModelAlphabetProblem::Alphabet(requirement));
  }
  model
    .frequencies()
    .filter(|pi| pi.len() != n_states)
    .map(|pi| ModelAlphabetProblem::Frequencies {
      given: pi.len(),
      n_states,
    })
}

enum ModelAlphabetProblem {
  Alphabet(AlphabetRequirement),
  Frequencies { given: usize, n_states: usize },
}

impl ModelAlphabetProblem {
  fn message(&self, partition: Option<&str>, model: &GtrModel, alphabet: AlphabetName) -> String {
    let model_name = model.name().to_string().to_lowercase();
    match self {
      Self::Alphabet(AlphabetRequirement {
        alphabet: required,
        reason,
      }) => {
        let subject = partition.map_or_else(
          || format!("the model {model_name} needs"),
          |label| format!("{label} has the model {model_name}, which needs"),
        );
        format!("{subject} the alphabet {required} because {reason}, but the alphabet is {alphabet}")
      },
      Self::Frequencies { given, n_states } => {
        let values = counted(*given, "value");
        let subject = partition.map_or_else(
          || format!("the model {model_name} has {values} in pi"),
          |label| format!("{label} has the model {model_name} with {values} in pi"),
        );
        format!(
          "{subject}, but the alphabet {alphabet} has {}, so pi needs {n_states}",
          counted(*n_states, "state")
        )
      },
    }
  }
}

fn layout_error<T>(message: String) -> Result<T, Report> {
  Err(Report::new(PartitionError {
    kind: PartitionErrorKind::Layout,
    message,
  }))
}

fn model_error<T>(message: String) -> Result<T, Report> {
  Err(Report::new(PartitionError {
    kind: PartitionErrorKind::Model,
    message,
  }))
}

#[derive(Clone, Debug, PartialEq, Eq, Display, Error)]
#[display("{message}")]
pub struct PartitionError {
  pub kind: PartitionErrorKind,
  pub message: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PartitionErrorKind {
  Layout,
  Model,
}
