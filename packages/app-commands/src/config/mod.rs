#[cfg(test)]
mod __tests__;

#[cfg(feature = "clap")]
pub mod catalog;
#[cfg(feature = "clap")]
pub mod cli_flags;
#[cfg(feature = "clap")]
pub mod cli_rules;
#[cfg(feature = "clap")]
pub mod code;
pub mod labels;
pub mod load;
pub mod properties;
pub mod resolve_paths;
pub mod schema;
pub mod schema_check;
pub mod settings;
pub mod source;
pub mod suggest;
