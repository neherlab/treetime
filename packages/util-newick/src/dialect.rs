use eyre::{Report, eyre};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NewickStructure {
  #[default]
  Classic,
  ENewick,
  Rich,
}

impl NewickStructure {
  pub const ALL: [Self; 3] = [Self::Classic, Self::ENewick, Self::Rich];

  pub const fn hybrid_tags(self) -> bool {
    matches!(self, Self::ENewick | Self::Rich)
  }

  pub const fn field_count(self) -> usize {
    match self {
      Self::Rich => 3,
      Self::Classic | Self::ENewick => 1,
    }
  }

  pub const fn weight(self) -> bool {
    matches!(self, Self::Rich)
  }

  pub const fn name(self) -> &'static str {
    match self {
      Self::Classic => "classic",
      Self::ENewick => "enewick",
      Self::Rich => "rich",
    }
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NewickAnnotations {
  #[default]
  Plain,
  Beast,
  Nhx,
  MrBayes,
}

impl NewickAnnotations {
  pub const ALL: [Self; 4] = [Self::Plain, Self::Beast, Self::Nhx, Self::MrBayes];

  pub const fn reserves_annotations(self) -> bool {
    !matches!(self, Self::Plain)
  }

  pub const fn name(self) -> &'static str {
    match self {
      Self::Plain => "plain",
      Self::Beast => "beast",
      Self::Nhx => "nhx",
      Self::MrBayes => "mrbayes",
    }
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NewickDialect {
  pub structure: NewickStructure,
  pub annotations: NewickAnnotations,
}

impl NewickDialect {
  pub const CLASSIC: Self = Self::new(NewickStructure::Classic, NewickAnnotations::Plain);
  pub const BEAST: Self = Self::new(NewickStructure::Classic, NewickAnnotations::Beast);
  pub const NHX: Self = Self::new(NewickStructure::Classic, NewickAnnotations::Nhx);
  pub const MRBAYES: Self = Self::new(NewickStructure::Classic, NewickAnnotations::MrBayes);
  pub const ENEWICK: Self = Self::new(NewickStructure::ENewick, NewickAnnotations::Plain);
  pub const RICH: Self = Self::new(NewickStructure::Rich, NewickAnnotations::Plain);
  pub const ENEWICK_BEAST: Self = Self::new(NewickStructure::ENewick, NewickAnnotations::Beast);

  pub const fn new(structure: NewickStructure, annotations: NewickAnnotations) -> Self {
    Self { structure, annotations }
  }

  pub fn pairs() -> impl Iterator<Item = Self> {
    NewickStructure::ALL.into_iter().flat_map(|structure| {
      NewickAnnotations::ALL
        .into_iter()
        .map(move |annotations| Self::new(structure, annotations))
    })
  }

  pub const fn rooting(self) -> bool {
    matches!(self.annotations, NewickAnnotations::Beast | NewickAnnotations::MrBayes)
      || matches!(self.structure, NewickStructure::Rich)
  }
}

impl fmt::Display for NewickDialect {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "{},{}", self.structure.name(), self.annotations.name())
  }
}

impl FromStr for NewickDialect {
  type Err = Report;

  fn from_str(text: &str) -> Result<Self, Self::Err> {
    let structure_names = NewickStructure::ALL.map(NewickStructure::name).join(", ");
    let annotation_names = NewickAnnotations::ALL.map(NewickAnnotations::name).join(", ");
    let error = || {
      eyre!(
        "{text:?} is not a Newick dialect: expected <structure>,<annotations> with structure one of {structure_names} and annotations one of {annotation_names}"
      )
    };
    let (structure, annotations) = text.split_once(',').ok_or_else(error)?;
    let structure = NewickStructure::ALL
      .into_iter()
      .find(|candidate| candidate.name() == structure)
      .ok_or_else(error)?;
    let annotations = NewickAnnotations::ALL
      .into_iter()
      .find(|candidate| candidate.name() == annotations)
      .ok_or_else(error)?;
    Ok(Self::new(structure, annotations))
  }
}
