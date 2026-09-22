//! The record each workspace crate leaves for the `pub-unused-report`
//! post-pass: the `pub` items it defines and the workspace items it uses.
//! Shared with the report binary by path, so it depends on serde only.

use std::collections::BTreeSet;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

/// Directory the lint writes records into and the report reads them from.
/// Collection is off when unset.
pub const RECORDS_DIR_ENV: &str = "TREETIME_LINTS_PUB_UNUSED_DIR";

/// Extension of a record file.
pub const RECORD_EXTENSION: &str = "json";

/// One compiled crate: a lib, bin, example, bench, or build script target.
#[derive(Serialize, Deserialize)]
pub struct CrateRecord {
    /// Absolute path of the crate root file, the target's `src_path` in
    /// `cargo metadata`.
    pub src_path: PathBuf,
    pub crate_name: String,
    /// A binary, example, or build script: nothing else can use its items.
    pub executable: bool,
    /// Reachable `pub` items the crate defines.
    pub defs: Vec<Def>,
    /// Keys of the crate's own items it uses.
    pub local_refs: BTreeSet<String>,
    /// Keys of the other workspace libraries' items it uses.
    pub foreign_refs: BTreeSet<String>,
}

/// One reachable `pub` item.
#[derive(Clone, Serialize, Deserialize)]
pub struct Def {
    /// Crate name plus definition path, the same from every crate:
    /// `treetime_io::csv::read`, `treetime::{impl#3}::eql`.
    pub key: String,
    /// Absolute path of the file that holds the item's name.
    pub file: PathBuf,
    /// The item's name as byte offsets into `file`.
    pub lo: u32,
    pub hi: u32,
    /// "function", "struct", ...
    pub descr: String,
    /// `crate::module::Item`, for the message.
    pub path: String,
    /// Key of the trait or type an associated item belongs to: when that is
    /// reported, the item is not.
    pub parent: Option<String>,
}
