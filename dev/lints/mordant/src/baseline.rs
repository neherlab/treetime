//! Ratchet mode. A baseline file records the accepted finding count per
//! (lint, file); runs suppress up to that many findings and surface only the
//! overflow, so mordant can gate CI on a brownfield codebase from day one.
//!
//! Regeneration: `MORDANT_BASELINE_WRITE=1 cargo dylint --all` emits nothing
//! and rewrites each compiled crate's section instead. Sections are keyed by
//! crate so parallel rustc processes never clobber another crate's entries;
//! the file itself is serialized with an exclusive file lock.
//!
//! Counts, not spans: line numbers drift with every edit, so a per-file count
//! is the only key that survives normal development. Moving a finding between
//! files consumes allowance in one file and overflows in the other, which is
//! the desired ratchet behavior.
//!
//! Severity: with a baseline configured, the baseline decides what fails the
//! run, not the lint level. A finding over the recorded count is printed as a
//! plain warning through the session's diagnostic context rather than as a
//! lint, so `-D warnings`, `[lints] warnings = "deny"` and `--cap-lints`
//! cannot turn it into an error that stops the crate and hides every finding
//! after it. `#[allow]` and `#[expect]` are still read, at the same node the
//! lint path reads them. Each crate that goes over prints one summary line and
//! appends itself to `target/mordant/over-baseline.txt`, which is the file CI
//! tests. Without a baseline nothing here applies and findings are ordinary
//! lints at their ordinary levels.
//!
//! Every diagnostic a mordant lint produces goes through one of the three
//! entry points here (`emit`, `emit_with_note`, `emit_hir_then`), so each
//! finding is weighed against the baseline exactly once.

use std::collections::BTreeMap;
use std::collections::HashMap;
use std::io::{Read, Seek, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};

use clippy_utils::diagnostics::{span_lint_and_help, span_lint_and_then, span_lint_hir_and_then};
use rustc_errors::Diag;
use rustc_hir::HirId;
use rustc_lint::{LateContext, LateLintPass, Lint, LintContext};
use rustc_span::{FileName, Span};

/// (lint, file relative to the workspace root): the unit the baseline counts.
type Key = (String, String);

/// What this process does with each finding; fixed for the whole run at
/// `setup`, since the environment variable that selects it cannot change
/// while rustc is running.
enum Mode {
    /// Stay silent up to the recorded count for each key and report the
    /// overflow as warnings.
    Ratchet {
        recorded: HashMap<Key, usize>,
        /// Findings weighed so far this run, per key.
        seen: Mutex<HashMap<Key, usize>>,
        /// Over-baseline findings reported (not allowed away) in this crate.
        over: AtomicUsize,
        /// Where a crate that went over appends its name and count.
        status_file: PathBuf,
    },
    /// Emit nothing; collect every finding and rewrite this crate's section of
    /// the file at `path` from `BaselineWriter::check_crate_post`.
    Record {
        path: PathBuf,
        recorded: Mutex<Vec<Key>>,
    },
}

pub struct Baseline {
    root: PathBuf,
    mode: Mode,
}

static STATE: OnceLock<Option<Baseline>> = OnceLock::new();

fn state() -> Option<&'static Baseline> {
    STATE.get().and_then(Option::as_ref)
}

fn write_mode() -> bool {
    std::env::var_os("MORDANT_BASELINE_WRITE").is_some()
}

type Doc = BTreeMap<String, BTreeMap<String, u64>>;

fn read_doc(bytes: &str) -> Doc {
    toml::from_str(bytes).unwrap_or_default()
}

/// Called once from `register_lints`.
pub fn setup(file_name: &Option<String>) {
    let _ = STATE.set(file_name.as_ref().and_then(|f| init(f)));
}

fn init(file_name: &str) -> Option<Baseline> {
    let record = write_mode();
    let mut dir = PathBuf::from(std::env::var_os("CARGO_MANIFEST_DIR")?);
    loop {
        let cand = dir.join(file_name);
        // In write mode the file may not exist yet; anchor at the workspace
        // root, which is where dylint.toml (the config that named us) lives.
        if cand.exists() || (record && dir.join("dylint.toml").exists()) {
            let mode = if record {
                Mode::Record {
                    path: cand,
                    recorded: Mutex::new(Vec::new()),
                }
            } else {
                Mode::Ratchet {
                    recorded: read_recorded(&cand),
                    seen: Mutex::new(HashMap::new()),
                    over: AtomicUsize::new(0),
                    status_file: status_file(&dir, cargo_target_dir().as_deref()),
                }
            };
            return Some(Baseline { root: dir, mode });
        }
        if !dir.pop() {
            return None;
        }
    }
}

/// `CARGO_TARGET_DIR` when it names a directory; cargo treats an empty
/// value as unset, so that is `None` here too.
fn cargo_target_dir() -> Option<PathBuf> {
    std::env::var_os("CARGO_TARGET_DIR")
        .filter(|d| !d.is_empty())
        .map(PathBuf::from)
}

/// `${CARGO_TARGET_DIR or <root>/target}`, a relative `CARGO_TARGET_DIR`
/// taken from the workspace root as cargo does.
pub(crate) fn target_dir(root: &Path) -> PathBuf {
    resolve_target_dir(root, cargo_target_dir().as_deref())
}

fn resolve_target_dir(root: &Path, target_dir: Option<&Path>) -> PathBuf {
    match target_dir {
        Some(dir) => root.join(dir),
        None => root.join("target"),
    }
}

/// `<target>/mordant/over-baseline.txt`.
fn status_file(root: &Path, target_dir: Option<&Path>) -> PathBuf {
    resolve_target_dir(root, target_dir)
        .join("mordant")
        .join("over-baseline.txt")
}

/// Sums every crate section of the file, since a (lint, file) key can appear
/// under more than one crate (a file shared by a lib and a bin target).
fn read_recorded(path: &Path) -> HashMap<Key, usize> {
    let doc = std::fs::read_to_string(path)
        .map(|s| read_doc(&s))
        .unwrap_or_default();
    let mut counts: HashMap<Key, usize> = HashMap::new();
    for section in doc.values() {
        for (key, n) in section {
            if let Some((lint, file)) = key.split_once(':') {
                *counts
                    .entry((lint.to_string(), file.to_string()))
                    .or_default() += *n as usize;
            }
        }
    }
    counts
}

fn rel_file(cx: &LateContext<'_>, b: &Baseline, span: Span) -> Option<String> {
    let FileName::Real(real) = cx.tcx.sess.source_map().span_to_filename(span) else {
        return None;
    };
    let path = real.local_path()?.to_path_buf();
    let rel = path.strip_prefix(&b.root).unwrap_or(&path);
    Some(rel.to_string_lossy().into_owned())
}

/// What the baseline decided about one finding.
enum Verdict {
    /// No baseline governs it: an ordinary lint at its ordinary level.
    Lint,
    /// Allowed or expected at its node, recorded in write mode, or within
    /// the recorded count: nothing is printed.
    Silent,
    /// Over the recorded count: printed as a warning no lint level can raise.
    Over(Over),
}

struct Over {
    lint: &'static Lint,
    recorded: usize,
    file: String,
    count: &'static AtomicUsize,
}

impl Over {
    fn report(
        self,
        cx: &LateContext<'_>,
        span: Span,
        msg: String,
        decorate: impl FnOnce(&mut Diag<'_, ()>),
    ) {
        self.count.fetch_add(1, Ordering::Relaxed);
        let mut diag = cx.tcx.dcx().struct_span_warn(span, msg);
        decorate(&mut diag);
        diag.note(format!(
            "`{}` over the mordant baseline ({} recorded for {})",
            self.lint.name_lower(),
            self.recorded,
            self.file,
        ));
        diag.emit();
    }
}

/// Weighs one finding against the baseline; each finding must pass through
/// here exactly once, since it counts. `hir_id` is the node the lint path
/// would read the level at, so the baseline honours exactly the `#[allow]` /
/// `#[expect]` the lint would have: a finding allowed there prints nothing
/// without a baseline, so it is neither recorded nor counted with one, and
/// an expectation there is fulfilled whatever the count says, as the lint did
/// fire under it.
fn weigh(cx: &LateContext<'_>, lint: &'static Lint, span: Span, hir_id: HirId) -> Verdict {
    let Some(b) = state() else {
        return Verdict::Lint;
    };
    let Some(file) = rel_file(cx, b, span) else {
        return Verdict::Lint;
    };
    let spec = cx.tcx.lint_level_spec_at_node(lint, hir_id);
    if let Some(expectation) = spec.lint_id() {
        cx.fulfill_expectation(expectation);
    }
    if spec.is_allow() || spec.is_expect() || span.in_external_macro(cx.tcx.sess.source_map()) {
        return Verdict::Silent;
    }
    let key = (lint.name_lower(), file);
    match &b.mode {
        Mode::Record { recorded, .. } => {
            recorded.lock().unwrap().push(key);
            Verdict::Silent
        }
        Mode::Ratchet {
            recorded,
            seen,
            over,
            ..
        } => {
            let limit = recorded.get(&key).copied().unwrap_or(0);
            let within = {
                let mut seen = seen.lock().unwrap();
                let n = seen.entry(key.clone()).or_default();
                *n += 1;
                *n <= limit
            };
            if within {
                return Verdict::Silent;
            }
            Verdict::Over(Over {
                lint,
                recorded: limit,
                file: key.1,
                count: over,
            })
        }
    }
}

/// `a`, `a and b`, `a, b and c`, with `conj` as the last separator.
pub fn join(items: &[String], conj: &str) -> String {
    match items {
        [] => String::new(),
        [one] => one.clone(),
        [head @ .., last] => format!("{} {conj} {last}", head.join(", ")),
    }
}

/// Every mordant lint reports through here so the ratchet sees all of them.
pub fn emit(
    cx: &LateContext<'_>,
    lint: &'static Lint,
    span: Span,
    msg: impl Into<String>,
    help: impl Into<String>,
) {
    let help: String = help.into();
    match weigh(cx, lint, span, cx.last_node_with_lint_attrs) {
        Verdict::Silent => {}
        Verdict::Lint => span_lint_and_help(cx, lint, span, msg.into(), None, help),
        Verdict::Over(over) => over.report(cx, span, msg.into(), |diag| {
            diag.help(help);
        }),
    }
}

/// `emit` with a secondary span: the finding is at `span`, the evidence for
/// it (the check, the other acquisition) is at `note_span`.
pub fn emit_with_note(
    cx: &LateContext<'_>,
    lint: &'static Lint,
    span: Span,
    msg: impl Into<String>,
    note_span: Span,
    note: impl Into<String>,
    help: impl Into<String>,
) {
    let (note, help): (String, String) = (note.into(), help.into());
    let decorate = |diag: &mut Diag<'_, ()>| {
        diag.span_note(note_span, note);
        diag.help(help);
    };
    match weigh(cx, lint, span, cx.last_node_with_lint_attrs) {
        Verdict::Silent => {}
        Verdict::Lint => span_lint_and_then(cx, lint, span, msg.into(), decorate),
        Verdict::Over(over) => over.report(cx, span, msg.into(), decorate),
    }
}

/// `emit` for a finding whose lint level is read at `hir_id` rather than at
/// the enclosing item, so an `#[allow]` on that node alone (a match arm, say)
/// is honored; `decorate` attaches the suggestion or help.
pub fn emit_hir_then(
    cx: &LateContext<'_>,
    lint: &'static Lint,
    hir_id: HirId,
    span: Span,
    msg: impl Into<String>,
    decorate: impl FnOnce(&mut Diag<'_, ()>),
) {
    match weigh(cx, lint, span, hir_id) {
        Verdict::Silent => {}
        Verdict::Lint => span_lint_hir_and_then(cx, lint, hir_id, span, msg.into(), decorate),
        Verdict::Over(over) => over.report(cx, span, msg.into(), decorate),
    }
}

/// The baseline section name for this compilation: the crate, with the bin
/// target appended, since one crate name can cover a lib and several bins.
fn section_name(cx: &LateContext<'_>) -> String {
    let name = cx
        .tcx
        .crate_name(rustc_hir::def_id::LOCAL_CRATE)
        .to_string();
    match std::env::var_os("CARGO_BIN_NAME") {
        Some(bin) => format!("{name} (bin {})", bin.to_string_lossy()),
        None => name,
    }
}

// Registered last: flushes write-mode recordings for this crate into the
// baseline file, replacing only this crate's section, or in ratchet mode
// prints the crate's over-baseline summary. It declares no lint of its own;
// rustc always runs a lintless pass, so nothing can `allow` it away.
rustc_session::declare_lint_pass!(BaselineWriter => []);

impl<'tcx> LateLintPass<'tcx> for BaselineWriter {
    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        match state() {
            None => {}
            Some(Baseline {
                mode: Mode::Ratchet {
                    over, status_file, ..
                },
                ..
            }) => summarize(cx, over.load(Ordering::Relaxed), status_file),
            Some(Baseline {
                mode: Mode::Record { path, recorded },
                ..
            }) => write_section(cx, path, recorded),
        }
    }
}

/// One line for the reader and one for CI. The status file is appended to,
/// never truncated: every crate is its own rustc process, so no process knows
/// it is the first. CI removes the file before the run and tests it is empty
/// or absent after.
fn summarize(cx: &LateContext<'_>, over: usize, status_file: &Path) {
    if over == 0 {
        return;
    }
    let name = section_name(cx);
    cx.tcx.dcx().warn(format!(
        "mordant: {over} finding(s) over the baseline in {name}"
    ));
    if let Some(dir) = status_file.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(status_file)
    else {
        return;
    };
    let _ = f.lock();
    let _ = f.write_all(format!("{name} {over}\n").as_bytes());
    let _ = f.unlock();
}

fn write_section(cx: &LateContext<'_>, path: &Path, recorded: &Mutex<Vec<Key>>) {
    let recorded: Vec<Key> = std::mem::take(&mut *recorded.lock().unwrap());
    let mut section: BTreeMap<String, u64> = BTreeMap::new();
    for (lint, file) in recorded {
        *section.entry(format!("{lint}:{file}")).or_default() += 1;
    }
    let name = section_name(cx);
    let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        // Not `truncate`: the file is read first, then rewritten in place
        // under the lock via `set_len(0)`.
        .truncate(false)
        .read(true)
        .write(true)
        .open(path)
    else {
        return;
    };
    // Parallel rustc processes write concurrently; the lock serializes the
    // read-modify-write.
    let _ = f.lock();
    let mut existing = String::new();
    let _ = f.read_to_string(&mut existing);
    let mut doc = read_doc(&existing);
    if section.is_empty() {
        doc.remove(&name);
    } else {
        doc.insert(name, section);
    }
    if let Ok(out) = toml::to_string_pretty(&doc) {
        let _ = f.set_len(0);
        let _ = f.rewind();
        let _ = f.write_all(out.as_bytes());
    }
    let _ = f.unlock();
}

#[cfg(test)]
mod tests {
    use super::status_file;
    use std::path::{Path, PathBuf};

    #[test]
    fn status_file_defaults_to_target_under_the_workspace_root() {
        assert_eq!(
            status_file(Path::new("/ws"), None),
            PathBuf::from("/ws/target/mordant/over-baseline.txt"),
        );
    }

    #[test]
    fn status_file_follows_cargo_target_dir() {
        assert_eq!(
            status_file(Path::new("/ws"), Some(Path::new("/elsewhere/tgt"))),
            PathBuf::from("/elsewhere/tgt/mordant/over-baseline.txt"),
        );
        assert_eq!(
            status_file(Path::new("/ws"), Some(Path::new("build/cargo"))),
            PathBuf::from("/ws/build/cargo/mordant/over-baseline.txt"),
        );
    }
}
