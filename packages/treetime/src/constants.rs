pub const TINY_NUMBER: f64 = 1e-12;

/// Floor added to transition matrices to prevent zero joint probabilities in mutation counting
pub const SUPERTINY_NUMBER: f64 = 1e-24;

/// fraction of length 'one_mutation' that is used as lower cut-off for branch lengths in GTR
pub const MIN_BRANCH_LENGTH_FRACTION: f64 = 1e-3;
