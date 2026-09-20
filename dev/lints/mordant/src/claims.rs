//! Shared machinery for lints that treat prose as checkable claims
//! (`stale_safety_comment`, `stale_panic_message`). A claim is a backticked
//! name in a comment or message; it holds if the name occurs in the file's
//! code with comments removed, or is defined by this crate or any crate it
//! links. Both callers run the same check: [`backticked_idents`] over the
//! text, then [`word_in`] against [`file_code_only`] and [`DefNames::contains`]
//! for each name.

use std::cell::OnceCell;
use std::collections::HashSet;

use rustc_hir::def_id::{DefId, DefIndex};
use rustc_lint::LateContext;
use rustc_metadata::creader::CStore;
use rustc_span::{BytePos, Span};

/// Keywords, primitives, and prelude names. Every one of them is defined
/// everywhere, so a claim naming one can never go stale and is not extracted.
pub(crate) const IDENT_STOPLIST: &[&str] = &[
    "self", "Self", "mut", "true", "false", "None", "Some", "Ok", "Err", "Vec", "Box", "drop",
    "u8", "u16", "u32", "u64", "u128", "usize", "i8", "i16", "i32", "i64", "i128", "isize", "f32",
    "f64", "bool", "str", "String", "unsafe", "SAFETY",
];

/// A backticked `name.ext` with one of these extensions is a file, not a
/// path expression whose last segment is an identifier.
const FILE_EXTENSIONS: &[&str] = &[
    "rs", "toml", "md", "txt", "json", "yaml", "yml", "lock", "c", "h", "cc", "cpp", "hpp", "zig",
    "js", "ts", "py", "go",
];

/// Backticked mentions in `text`, reduced to their final identifier:
/// `` `self.frames` `` and `` `VM::wake()` `` yield `frames` and `wake`.
/// Skipped: file paths and names with a [`FILE_EXTENSIONS`] suffix, anything
/// not identifier-shaped after the reduction, and the [`IDENT_STOPLIST`].
/// Duplicates are kept.
pub(crate) fn backticked_idents(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(start) = rest.find('`') {
        let after = &rest[start + 1..];
        let Some(end) = after.find('`') else { break };
        let token = &after[..end];
        rest = &after[end + 1..];
        if token.contains('/')
            || token
                .rsplit_once('.')
                .is_some_and(|(_, ext)| FILE_EXTENSIONS.contains(&ext))
        {
            continue;
        }
        let last = token
            .trim_end_matches("()")
            .rsplit(&[':', '.'][..])
            .next()
            .unwrap_or(token);
        let is_ident = !last.is_empty()
            && last
                .chars()
                .next()
                .is_some_and(|c| c.is_alphabetic() || c == '_')
            && last.chars().all(|c| c.is_alphanumeric() || c == '_');
        if is_ident && !IDENT_STOPLIST.contains(&last) {
            out.push(last.to_string());
        }
    }
    out
}

/// The contiguous run of `//` lines directly above `span`'s first line (a
/// blank line or code ends the run; nothing directly above means `None`),
/// with the slashes stripped and the lines joined by `\n`, plus the span
/// covering exactly those lines, which is what the callers report on.
pub(crate) fn comment_above(cx: &LateContext<'_>, span: Span) -> Option<(String, Span)> {
    let sm = cx.tcx.sess.source_map();
    let lines = sm.span_to_lines(span).ok()?;
    let file = lines.file;
    let first = lines.lines.first()?.line_index;
    let mut collected: Vec<String> = Vec::new();
    let mut top = first;
    for i in (0..first).rev() {
        let line = file.get_line(i)?;
        let trimmed = line.trim_start();
        if trimmed.starts_with("//") {
            collected.push(trimmed.trim_start_matches('/').to_string());
            top = i;
        } else {
            break;
        }
    }
    if collected.is_empty() {
        return None;
    }
    collected.reverse();
    let lo = file.line_bounds(top).start;
    let hi = file.line_bounds(first.saturating_sub(1)).end;
    Some((
        collected.join("\n"),
        Span::with_root_ctxt(BytePos(lo.0), BytePos(hi.0)),
    ))
}

/// The source of `span`'s whole file with the rest of every line cut at its
/// first `//`, so a comment can never vouch for its own claims. The cut is
/// textual: a `//` inside a string literal loses the rest of that line too,
/// and `/* */` comments are left in.
pub(crate) fn file_code_only(cx: &LateContext<'_>, span: Span) -> String {
    let file_src = cx
        .tcx
        .sess
        .source_map()
        .span_to_lines(span)
        .ok()
        .and_then(|l| l.file.src.as_ref().map(|s| s.to_string()))
        .unwrap_or_default();
    file_src
        .lines()
        .map(|l| l.split("//").next().unwrap_or(""))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Whether `ident` occurs in `haystack` as a whole identifier token, so
/// `frame` does not match `frames` or `frame_count`.
pub(crate) fn word_in(haystack: &str, ident: &str) -> bool {
    haystack
        .split(|c: char| !(c.is_alphanumeric() || c == '_'))
        .any(|w| w == ident)
}

/// The definition names a claim may point at beyond its own file: everything
/// this crate defines, and everything any linked crate defines. Built once
/// per crate in the callers' `check_crate`. The upstream set is large (std
/// alone is tens of thousands of names) and a crate whose claims all resolve
/// locally never needs it, so it is built on the first lookup that misses
/// locally.
pub(crate) struct DefNames {
    local: HashSet<String>,
    upstream: OnceCell<HashSet<String>>,
}

impl DefNames {
    pub(crate) fn collect(cx: &LateContext<'_>) -> Self {
        Self {
            local: crate_def_names(cx),
            upstream: OnceCell::new(),
        }
    }

    pub(crate) fn contains(&self, cx: &LateContext<'_>, ident: &str) -> bool {
        self.local.contains(ident)
            || self
                .upstream
                .get_or_init(|| upstream_def_names(cx))
                .contains(ident)
    }
}

/// Every named definition in every crate this one links, std included:
/// items, variants, and fields alike, since all have a `def_key` name. A
/// SAFETY comment in a workspace member routinely names the sibling crate or
/// std type that owns the invariant, and none of those are in the local HIR.
fn upstream_def_names(cx: &LateContext<'_>) -> HashSet<String> {
    let mut names = HashSet::new();
    let cstore = CStore::from_tcx(cx.tcx);
    for &cnum in cx.tcx.crates(()) {
        for i in 0..cstore.num_def_ids_untracked(cnum) {
            let did = DefId {
                krate: cnum,
                index: DefIndex::from_usize(i),
            };
            // A proc-macro crate's def table is sparse: only the root and the
            // macros themselves are encoded, and `def_key` unwraps on the
            // holes. An absent entry decodes its path hash as zero, which no
            // real definition has, so that is the probe for a hole.
            if cx.tcx.def_path_hash(did).local_hash().as_u64() == 0 {
                continue;
            }
            if let Some(name) = cx.tcx.def_key(did).get_opt_name() {
                names.insert(name.to_string());
            }
        }
    }
    names
}

/// Every name this crate defines: each item, and every variant and field of
/// each local ADT, since a claim usually names the field it guards.
fn crate_def_names(cx: &LateContext<'_>) -> HashSet<String> {
    let mut names = HashSet::new();
    for def_id in cx.tcx.hir_crate_items(()).definitions() {
        if let Some(name) = cx.tcx.opt_item_name(def_id.to_def_id()) {
            names.insert(name.to_string());
        }
        if let rustc_hir::def::DefKind::Struct
        | rustc_hir::def::DefKind::Enum
        | rustc_hir::def::DefKind::Union = cx.tcx.def_kind(def_id)
        {
            let adt = cx.tcx.adt_def(def_id);
            for v in adt.variants() {
                names.insert(v.name.to_string());
                names.extend(v.fields.iter().map(|f| f.name.to_string()));
            }
        }
    }
    names
}
