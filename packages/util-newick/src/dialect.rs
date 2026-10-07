use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NewickDialect {
  Classic,
  Beast,
  MrBayes,
  Nhx,
  ENewick,
  Rich,
}

impl NewickDialect {
  pub const ALL: [Self; 6] = [
    Self::Rich,
    Self::ENewick,
    Self::Nhx,
    Self::MrBayes,
    Self::Beast,
    Self::Classic,
  ];

  pub const fn features(self) -> DialectFeatures {
    match self {
      Self::Classic => DialectFeatures {
        comments: CommentKind::Plain,
        reserves_annotations: false,
        hybrid_tags: false,
        rich_fields: false,
        rooting: false,
        weight: false,
      },
      Self::Beast => DialectFeatures {
        comments: CommentKind::Beast,
        reserves_annotations: true,
        hybrid_tags: false,
        rich_fields: false,
        rooting: true,
        weight: false,
      },
      Self::MrBayes => DialectFeatures {
        comments: CommentKind::MrBayes,
        reserves_annotations: true,
        hybrid_tags: false,
        rich_fields: false,
        rooting: true,
        weight: false,
      },
      Self::Nhx => DialectFeatures {
        comments: CommentKind::Nhx,
        reserves_annotations: true,
        hybrid_tags: false,
        rich_fields: false,
        rooting: false,
        weight: false,
      },
      Self::ENewick => DialectFeatures {
        comments: CommentKind::Plain,
        reserves_annotations: true,
        hybrid_tags: true,
        rich_fields: false,
        rooting: false,
        weight: false,
      },
      Self::Rich => DialectFeatures {
        comments: CommentKind::Plain,
        reserves_annotations: true,
        hybrid_tags: true,
        rich_fields: true,
        rooting: true,
        weight: true,
      },
    }
  }

  pub const fn name(self) -> &'static str {
    match self {
      Self::Classic => "classic",
      Self::Beast => "beast",
      Self::MrBayes => "mrbayes",
      Self::Nhx => "nhx",
      Self::ENewick => "enewick",
      Self::Rich => "rich",
    }
  }
}

impl fmt::Display for NewickDialect {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    f.write_str(self.name())
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DialectFeatures {
  pub comments: CommentKind,
  pub reserves_annotations: bool,
  pub hybrid_tags: bool,
  pub rich_fields: bool,
  pub rooting: bool,
  pub weight: bool,
}

impl DialectFeatures {
  pub const fn field_count(self) -> usize {
    if self.rich_fields { 3 } else { 1 }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CommentKind {
  Plain,
  Beast,
  Nhx,
  MrBayes,
}
