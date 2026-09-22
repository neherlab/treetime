pub mod cost_function;
pub(crate) mod find_best_root;
pub mod find_best_split;
pub(crate) mod method_brent;
pub(crate) mod method_golden_section;
pub(crate) mod method_grid_search;
pub mod params;

#[cfg(test)]
mod __tests__;
