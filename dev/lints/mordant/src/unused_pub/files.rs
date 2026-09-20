//! The two files each member writes under `<target>/mordant/unused_pub/`.
//! `<crate>.defs` has one line per `pub` item a library defines.
//! `<package>.<crate>.<lib|bin>.refs` has one line per workspace item the
//! crate uses. Items are keyed by crate name plus definition path, which
//! reads the same from every crate. A file is written whole under a
//! temporary name and renamed, so a reader never sees half of one.

use std::collections::{BTreeSet, HashSet};
use std::path::{Path, PathBuf};

use rustc_hir::def_id::DefId;

use super::workspace::Kind;
use rustc_lint::LateContext;
use rustc_middle::ty::TyCtxt;
use rustc_span::{BytePos, FileName, Span, SyntaxContext};

/// One `pub` item, as its defining crate describes it for the crate that
/// will print it.
pub struct Def {
    pub key: String,
    /// Source file as rustc named it, relative to the workspace root for
    /// a member.
    pub file: String,
    /// The item's name, as byte offsets into `file`.
    pub lo: u32,
    pub hi: u32,
    /// "function", "struct", ...
    pub descr: String,
    /// `crate::module::Item`, for the message.
    pub path: String,
    /// Key of the trait or type the item belongs to, or empty: when that
    /// is reported, the item is not.
    pub parent: String,
}

/// Crate name plus definition path: `bun_core::fmt::raw`,
/// `bun_css::{impl#3}::eql`. The same string whichever crate computes it.
pub fn key(tcx: TyCtxt<'_>, def_id: DefId) -> String {
    format!(
        "{}{}",
        tcx.crate_name(def_id.krate),
        tcx.def_path(def_id).to_string_no_crate_verbose()
    )
}

pub fn defs_path(dir: &Path, crate_name: &str) -> PathBuf {
    dir.join(format!("{crate_name}.defs"))
}

pub fn refs_path(dir: &Path, package: &str, crate_name: &str, kind: Kind) -> PathBuf {
    let kind = match kind {
        Kind::Lib => "lib",
        Kind::Bin => "bin",
    };
    dir.join(format!("{package}.{crate_name}.{kind}.refs"))
}

pub fn write_defs<'a>(path: &Path, defs: impl Iterator<Item = &'a Def>) {
    let mut out = String::new();
    for d in defs {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            d.key, d.file, d.lo, d.hi, d.descr, d.path, d.parent
        ));
    }
    write_whole(path, &out);
}

pub fn read_defs(path: &Path) -> Vec<Def> {
    let Ok(text) = std::fs::read_to_string(path) else {
        return Vec::new();
    };
    text.lines()
        .filter_map(|line| {
            let mut f = line.split('\t');
            Some(Def {
                key: f.next()?.into(),
                file: f.next()?.into(),
                lo: f.next()?.parse().ok()?,
                hi: f.next()?.parse().ok()?,
                descr: f.next()?.into(),
                path: f.next()?.into(),
                parent: f.next().unwrap_or("").into(),
            })
        })
        .collect()
}

pub fn write_refs(path: &Path, refs: &BTreeSet<String>) {
    let mut out = String::new();
    for r in refs {
        out.push_str(r);
        out.push('\n');
    }
    write_whole(path, &out);
}

/// Every key any `.defs` file in `dir` lists: the items a use is worth
/// recording for.
pub fn all_def_keys(dir: &Path) -> HashSet<String> {
    let mut keys = HashSet::new();
    for path in files_with_extension(dir, "defs") {
        keys.extend(read_defs(&path).into_iter().map(|d| d.key));
    }
    keys
}

/// The union of every `.refs` file in `dir`. A file whose package is no
/// longer a member is deleted instead, so a removed crate stops keeping
/// items in use.
pub fn all_refs(dir: &Path, members: &BTreeSet<String>) -> HashSet<String> {
    let mut refs = HashSet::new();
    for path in files_with_extension(dir, "refs") {
        let package = path
            .file_name()
            .and_then(|n| n.to_str())
            .and_then(|n| n.split('.').next())
            .unwrap_or("");
        if !members.contains(package) {
            let _ = std::fs::remove_file(&path);
            continue;
        }
        if let Ok(text) = std::fs::read_to_string(&path) {
            refs.extend(text.lines().map(str::to_string));
        }
    }
    refs
}

fn files_with_extension(dir: &Path, ext: &str) -> Vec<PathBuf> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    entries
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().and_then(|e| e.to_str()) == Some(ext))
        .collect()
}

fn write_whole(path: &Path, contents: &str) {
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let tmp = path.with_extension(format!("tmp{}", std::process::id()));
    if std::fs::write(&tmp, contents).is_ok() {
        let _ = std::fs::rename(&tmp, path);
    }
}

/// Where `span` starts and ends in its file, if that is a real file.
pub fn locate(cx: &LateContext<'_>, span: Span) -> Option<(String, u32, u32)> {
    let sm = cx.tcx.sess.source_map();
    let file = sm.lookup_source_file(span.lo());
    let FileName::Real(real) = &file.name else {
        return None;
    };
    let path = real.local_path()?.to_string_lossy().into_owned();
    Some((
        path,
        (span.lo() - file.start_pos).0,
        (span.hi() - file.start_pos).0,
    ))
}

/// A span at `lo..hi` of `file`, loading the file into this session's
/// source map so the diagnostic can quote it. `root` is tried as a base
/// when the path alone does not open.
pub fn span_in(cx: &LateContext<'_>, root: &Path, file: &str, lo: u32, hi: u32) -> Option<Span> {
    let sm = cx.tcx.sess.source_map();
    let loaded = sm
        .load_file(Path::new(file))
        .or_else(|_| sm.load_file(&root.join(file)))
        .ok()?;
    if loaded
        .src
        .as_ref()
        .is_none_or(|src| hi as usize > src.len())
    {
        return None;
    }
    Some(Span::new(
        loaded.start_pos + BytePos(lo),
        loaded.start_pos + BytePos(hi),
        SyntaxContext::root(),
        None,
    ))
}
