pub const BIG_NUMBER: f64 = 1e10;
pub const TINY_NUMBER: f64 = 1e-12;
pub const SUPERTINY_NUMBER: f64 = 1e-24;
pub const MIN_LOG: f64 = -1e8; // minimal log value
pub const MIN_BRANCH_LENGTH: f64 = 1e-3; // fraction of length 'one_mutation' that is used as lower cut-off for branch lengths in GTR
pub const OVER_DISPERSION: usize = 10;

// distribution parameters
pub const BRANCH_GRID_SIZE_ROUGH: usize = 75;
pub const NODE_GRID_SIZE_ROUGH: usize = 60;
pub const N_INTEGRAL_ROUGH: usize = 60;

pub const BRANCH_GRID_SIZE: usize = 125;
pub const NODE_GRID_SIZE: usize = 100;
pub const N_INTEGRAL: usize = 100;

pub const BRANCH_GRID_SIZE_FINE: usize = 200;
pub const NODE_GRID_SIZE_FINE: usize = 180;
pub const N_INTEGRAL_FINE: usize = 150;

pub const BRANCH_GRID_SIZE_ULTRA: usize = 300;
pub const NODE_GRID_SIZE_ULTRA: usize = 400;
pub const N_INTEGRAL_ULTRA: usize = 250;

// distribution parameters for FFT
pub const FFT_FWHM_GRID_SIZE: usize = 150;

pub const MIN_INTEGRATION_PEAK: f64 = 0.001;

// clocktree parameters
pub const BRANCH_LEN_PENALTY: usize = 0;
pub const MAX_BRANCH_LENGTH: f64 = 4.0; // only relevant for branch length optimization and time trees - upper boundary of interpolator objects
pub const NINTEGRAL: usize = 300;
pub const REL_TOL_PRUNE: f64 = 0.01;
pub const REL_TOL_REFINE: f64 = 0.05;
pub const NIQD: usize = 3;

// treetime
// autocorrelated molecular clock coefficients
pub const MU_ALPHA: f64 = 1.0;
pub const MU_BETA: f64 = 1.0;
