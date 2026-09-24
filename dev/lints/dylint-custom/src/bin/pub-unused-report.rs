//! Reports the reachable `pub` items that no workspace crate uses, from the
//! records the `pub_unused_in_workspace` lint leaves while `cargo dylint`
//! compiles the workspace. Runs after cargo exits, when every record is final.
//!
//! Usage: `pub-unused-report <workspace-manifest-path> [--exclude-crate <crate>]...`,
//! with the records directory in `TREETIME_LINTS_PUB_UNUSED_DIR` as for the lint.
//! An excluded crate publishes its API for users outside the workspace or serves
//! only test targets: its items are not reported, and its uses of other crates'
//! items still count.
//!
//! Findings fail the report with exit status 1. A workspace target
//! whose record is missing or unreadable makes the report incomplete: it is
//! named and the exit status is 1.

#[path = "../lints/pub_unused_in_workspace/record.rs"]
mod record;

use std::collections::{BTreeSet, HashSet};
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use cargo_metadata::{Metadata, MetadataCommand, TargetKind};

use record::{CrateRecord, Def, RECORD_EXTENSION, RECORDS_DIR_ENV};

fn main() -> ExitCode {
    let (Some(args), Some(records_dir)) = (parse_args(std::env::args_os().skip(1)), std::env::var_os(RECORDS_DIR_ENV))
    else {
        eprintln!(
            "usage: {RECORDS_DIR_ENV}=<records-dir> pub-unused-report <workspace-manifest-path> [--exclude-crate <crate>]..."
        );
        return ExitCode::from(2);
    };
    match run(Path::new(&records_dir), &args.manifest_path, &args.excluded_crates) {
        Ok(Completeness::Complete) => ExitCode::SUCCESS,
        Ok(Completeness::Findings) => ExitCode::FAILURE,
        Ok(Completeness::Incomplete) => ExitCode::FAILURE,
        Err(err) => {
            eprintln!("error: pub-unused-report: {err}");
            ExitCode::FAILURE
        },
    }
}

fn run(records_dir: &Path, manifest_path: &Path, excluded_crates: &BTreeSet<String>) -> Result<Completeness, String> {
    let metadata = MetadataCommand::new()
        .manifest_path(manifest_path)
        .no_deps()
        .exec()
        .map_err(|err| format!("cargo metadata for {}: {err}", manifest_path.display()))?;
    let targets = workspace_targets(&metadata);
    let loaded = load_records(records_dir)?;
    for (path, err) in &loaded.unreadable {
        eprintln!("error: pub-unused-report: cannot read record {}: {err}", path.display());
    }
    let records = select_records(&targets, loaded.records);
    let missing = missing_targets(&targets, &records);
    for target in &missing {
        eprintln!(
            "error: pub-unused-report: no record for {} `{}` of package `{}` ({})",
            target.kind,
            target.name,
            target.package,
            target.src_path.display()
        );
    }
    if !missing.is_empty() || !loaded.unreadable.is_empty() {
        eprintln!(
            "error: pub-unused-report: the report is incomplete and nothing is judged. A target has no \
             record when its compilation failed, when a compiler cache (RUSTC_WRAPPER) returned it \
             without running the lint, or when {RECORDS_DIR_ENV} was not set for `cargo dylint`"
        );
        return Ok(Completeness::Incomplete);
    }
    let root = metadata.workspace_root.as_std_path();
    let unused = unused_defs(&records, excluded_crates);
    for def in &unused {
        eprint!("{}", render(def, root));
    }
    if !unused.is_empty() {
        eprintln!(
            "error: pub_unused_in_workspace: {} public item(s) that no workspace crate uses",
            unused.len()
        );
        return Ok(Completeness::Findings);
    }
    Ok(Completeness::Complete)
}

#[derive(Debug, PartialEq)]
struct Args {
    manifest_path: PathBuf,
    excluded_crates: BTreeSet<String>,
}

/// The manifest path, then any number of `--exclude-crate <crate>` pairs.
fn parse_args(mut args: impl Iterator<Item = OsString>) -> Option<Args> {
    let manifest_path = PathBuf::from(args.next()?);
    let mut excluded_crates = BTreeSet::new();
    while let Some(flag) = args.next() {
        if flag != "--exclude-crate" {
            return None;
        }
        excluded_crates.insert(args.next()?.into_string().ok()?);
    }
    Some(Args {
        manifest_path,
        excluded_crates,
    })
}

#[derive(Debug, PartialEq)]
enum Completeness {
    Complete,
    Findings,
    Incomplete,
}

/// A workspace target whose compilation writes a record.
#[derive(Debug, PartialEq)]
struct WorkspaceTarget {
    package: String,
    name: String,
    kind: String,
    src_path: PathBuf,
    /// A bench or build script: `cargo check --all-targets` may compile a
    /// bench as a test harness, which writes no record, so its record is used
    /// when present and not required.
    required: bool,
}

/// Every non-test target of every workspace member.
fn workspace_targets(metadata: &Metadata) -> Vec<WorkspaceTarget> {
    metadata
        .workspace_packages()
        .into_iter()
        .flat_map(|package| {
            package.targets.iter().filter_map(|target| {
                let required = target.kind.iter().any(|kind| {
                    matches!(
                        kind,
                        TargetKind::Lib
                            | TargetKind::RLib
                            | TargetKind::DyLib
                            | TargetKind::CDyLib
                            | TargetKind::StaticLib
                            | TargetKind::ProcMacro
                            | TargetKind::Bin
                            | TargetKind::Example
                    )
                });
                let optional = target
                    .kind
                    .iter()
                    .any(|kind| matches!(kind, TargetKind::Bench | TargetKind::CustomBuild));
                (required || optional).then(|| WorkspaceTarget {
                    package: package.name.to_string(),
                    name: target.name.clone(),
                    kind: target.kind.iter().map(ToString::to_string).collect::<Vec<_>>().join(","),
                    src_path: target.src_path.clone().into_std_path_buf(),
                    required,
                })
            })
        })
        .collect()
}

struct LoadedRecords {
    records: Vec<CrateRecord>,
    unreadable: Vec<(PathBuf, String)>,
}

fn load_records(dir: &Path) -> Result<LoadedRecords, String> {
    let entries = match std::fs::read_dir(dir) {
        Ok(entries) => entries,
        Err(err) if err.kind() == std::io::ErrorKind::NotFound => {
            return Ok(LoadedRecords {
                records: vec![],
                unreadable: vec![],
            });
        },
        Err(err) => return Err(format!("cannot list {}: {err}", dir.display())),
    };
    let mut loaded = LoadedRecords {
        records: vec![],
        unreadable: vec![],
    };
    for entry in entries {
        let path = entry.map_err(|err| format!("cannot list {}: {err}", dir.display()))?.path();
        if path.extension().is_none_or(|ext| ext != RECORD_EXTENSION) {
            continue;
        }
        let parsed = std::fs::read(&path)
            .map_err(|err| err.to_string())
            .and_then(|bytes| serde_json::from_slice::<CrateRecord>(&bytes).map_err(|err| err.to_string()));
        match parsed {
            Ok(record) => loaded.records.push(record),
            Err(err) => loaded.unreadable.push((path, err)),
        }
    }
    Ok(loaded)
}

/// The records of current workspace targets. Records of removed targets or
/// packages are left out, so they neither keep items in use nor report any.
fn select_records(targets: &[WorkspaceTarget], records: Vec<CrateRecord>) -> Vec<CrateRecord> {
    let current = targets.iter().map(|target| target.src_path.as_path()).collect::<HashSet<_>>();
    records
        .into_iter()
        .filter(|record| current.contains(record.src_path.as_path()))
        .collect()
}

fn missing_targets<'a>(targets: &'a [WorkspaceTarget], records: &[CrateRecord]) -> Vec<&'a WorkspaceTarget> {
    let recorded = records.iter().map(|record| record.src_path.as_path()).collect::<HashSet<_>>();
    targets
        .iter()
        .filter(|target| target.required && !recorded.contains(target.src_path.as_path()))
        .collect()
}

/// Items no record uses, in file order. A library item counts as used when
/// its own crate or any other crate uses it; an executable's item only when
/// its own crate does, since nothing can depend on an executable. An item
/// whose trait or type is reported itself is left out, and so is every item of
/// an excluded crate.
fn unused_defs<'a>(records: &'a [CrateRecord], excluded_crates: &BTreeSet<String>) -> Vec<&'a Def> {
    let foreign = records
        .iter()
        .flat_map(|record| record.foreign_refs.iter().map(String::as_str))
        .collect::<HashSet<_>>();
    let mut unused = records
        .iter()
        .filter(|record| !excluded_crates.contains(&record.crate_name))
        .flat_map(|record| {
            let foreign = &foreign;
            record.defs.iter().filter(move |def| {
                !record.local_refs.contains(&def.key) && (record.executable || !foreign.contains(def.key.as_str()))
            })
        })
        .collect::<Vec<_>>();
    let reported = unused.iter().map(|def| def.key.as_str()).collect::<BTreeSet<_>>();
    unused.retain(|def| def.parent.as_deref().is_none_or(|parent| !reported.contains(parent)));
    unused.sort_by(|a, b| a.file.cmp(&b.file).then(a.lo.cmp(&b.lo)));
    unused
}

fn render(def: &Def, root: &Path) -> String {
    let shown = def.file.strip_prefix(root).unwrap_or(&def.file);
    let position = std::fs::read(&def.file)
        .ok()
        .and_then(|src| line_col(&src, def.lo as usize))
        .map_or_else(String::new, |(line, col)| format!(":{line}:{col}"));
    format!(
        "warning: {} `{}` is public, but nothing in the workspace uses it\n  --> {}{position}\n   = help: remove it; if code under a `cfg`, target, or feature not compiled here uses it, gate the item the same way\n   = note: `#[warn(pub_unused_in_workspace)]` on by default\n\n",
        def.descr,
        def.path,
        shown.display(),
    )
}

/// 1-based line and column (in characters) of byte `offset` in `src`.
fn line_col(src: &[u8], offset: usize) -> Option<(usize, usize)> {
    let before = std::str::from_utf8(src.get(..offset)?).ok()?;
    let line = before.matches('\n').count() + 1;
    let line_start = before.rfind('\n').map_or(0, |pos| pos + 1);
    let col = before[line_start..].chars().count() + 1;
    Some((line, col))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pub_unused_report_library_item_used_by_other_crate_is_not_reported() {
        let records = vec![
            record("/ws/a/src/lib.rs", "a", false, &[def("a::f", "/ws/a/src/lib.rs", 10, None)], &[], &[]),
            record("/ws/b/src/main.rs", "b", true, &[], &[], &["a::f"]),
        ];
        let expected: Vec<&str> = vec![];
        assert_eq!(expected, keys(&unused_defs(&records, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_library_item_used_only_by_own_crate_is_not_reported() {
        let records = vec![record(
            "/ws/a/src/lib.rs",
            "a",
            false,
            &[def("a::f", "/ws/a/src/lib.rs", 10, None)],
            &["a::f"],
            &[],
        )];
        let expected: Vec<&str> = vec![];
        assert_eq!(expected, keys(&unused_defs(&records, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_library_item_used_by_nobody_is_reported() {
        let records = vec![
            record("/ws/a/src/lib.rs", "a", false, &[def("a::f", "/ws/a/src/lib.rs", 10, None)], &[], &[]),
            record("/ws/b/src/main.rs", "b", true, &[], &[], &["a::g"]),
        ];
        assert_eq!(vec!["a::f"], keys(&unused_defs(&records, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_executable_item_is_not_kept_by_same_key_in_another_crate() {
        let records = vec![
            record("/ws/cli/src/main.rs", "tool", true, &[def("tool::run", "/ws/cli/src/main.rs", 5, None)], &[], &[]),
            record("/ws/x/src/main.rs", "x", true, &[], &[], &["tool::run"]),
        ];
        assert_eq!(vec!["tool::run"], keys(&unused_defs(&records, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_member_of_reported_type_is_left_out() {
        let records = vec![record(
            "/ws/a/src/lib.rs",
            "a",
            false,
            &[
                def("a::S", "/ws/a/src/lib.rs", 10, None),
                def("a::{impl#0}::m", "/ws/a/src/lib.rs", 40, Some("a::S")),
            ],
            &[],
            &[],
        )];
        assert_eq!(vec!["a::S"], keys(&unused_defs(&records, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_member_of_used_type_is_reported() {
        let records = vec![
            record(
                "/ws/a/src/lib.rs",
                "a",
                false,
                &[
                    def("a::S", "/ws/a/src/lib.rs", 10, None),
                    def("a::{impl#0}::m", "/ws/a/src/lib.rs", 40, Some("a::S")),
                ],
                &[],
                &[],
            ),
            record("/ws/b/src/main.rs", "b", true, &[], &[], &["a::S"]),
        ];
        assert_eq!(vec!["a::{impl#0}::m"], keys(&unused_defs(&records, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_findings_are_in_file_order() {
        let records = vec![
            record("/ws/b/src/lib.rs", "b", false, &[def("b::g", "/ws/b/src/lib.rs", 5, None)], &[], &[]),
            record(
                "/ws/a/src/lib.rs",
                "a",
                false,
                &[def("a::y", "/ws/a/src/lib.rs", 50, None), def("a::x", "/ws/a/src/lib.rs", 7, None)],
                &[],
                &[],
            ),
        ];
        assert_eq!(vec!["a::x", "a::y", "b::g"], keys(&unused_defs(&records, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_records_of_removed_targets_are_left_out() {
        let targets = vec![target("/ws/a/src/lib.rs", true)];
        let records = vec![
            record("/ws/a/src/lib.rs", "a", false, &[def("a::f", "/ws/a/src/lib.rs", 10, None)], &[], &[]),
            record("/ws/gone/src/main.rs", "gone", true, &[], &[], &["a::f"]),
        ];
        let selected = select_records(&targets, records);
        assert_eq!(vec!["a::f"], keys(&unused_defs(&selected, &BTreeSet::new())));
    }

    #[test]
    fn pub_unused_report_missing_required_target_is_named() {
        let targets = vec![
            target("/ws/a/src/lib.rs", true),
            target("/ws/b/src/main.rs", true),
            target("/ws/a/benches/bench.rs", false),
        ];
        let records = vec![record("/ws/a/src/lib.rs", "a", false, &[], &[], &[])];
        let expected = vec![&targets[1]];
        assert_eq!(expected, missing_targets(&targets, &records));
    }

    #[test]
    fn pub_unused_report_excluded_crate_item_is_not_reported() {
        let records = vec![
            record("/ws/a/src/lib.rs", "a", false, &[def("a::f", "/ws/a/src/lib.rs", 10, None)], &[], &[]),
            record("/ws/b/src/lib.rs", "b", false, &[def("b::g", "/ws/b/src/lib.rs", 10, None)], &[], &[]),
        ];
        let excluded = BTreeSet::from(["a".to_owned()]);
        assert_eq!(vec!["b::g"], keys(&unused_defs(&records, &excluded)));
    }

    #[test]
    fn pub_unused_report_excluded_crate_use_keeps_item_used() {
        let records = vec![
            record("/ws/a/src/lib.rs", "a", false, &[], &[], &["b::g"]),
            record("/ws/b/src/lib.rs", "b", false, &[def("b::g", "/ws/b/src/lib.rs", 10, None)], &[], &[]),
        ];
        let excluded = BTreeSet::from(["a".to_owned()]);
        let expected: Vec<&str> = vec![];
        assert_eq!(expected, keys(&unused_defs(&records, &excluded)));
    }

    #[test]
    fn pub_unused_report_parse_args_reads_manifest_and_exclusions() {
        let expected = Some(Args {
            manifest_path: PathBuf::from("/ws/Cargo.toml"),
            excluded_crates: BTreeSet::from(["a".to_owned(), "b".to_owned()]),
        });
        assert_eq!(
            expected,
            parse_args(os_args(&["/ws/Cargo.toml", "--exclude-crate", "a", "--exclude-crate", "b"]))
        );
    }

    #[test]
    fn pub_unused_report_parse_args_reads_manifest_alone() {
        let expected = Some(Args {
            manifest_path: PathBuf::from("/ws/Cargo.toml"),
            excluded_crates: BTreeSet::new(),
        });
        assert_eq!(expected, parse_args(os_args(&["/ws/Cargo.toml"])));
    }

    #[test]
    fn pub_unused_report_parse_args_rejects_missing_manifest() {
        assert_eq!(None, parse_args(os_args(&[])));
    }

    #[test]
    fn pub_unused_report_parse_args_rejects_exclusion_without_crate() {
        assert_eq!(None, parse_args(os_args(&["/ws/Cargo.toml", "--exclude-crate"])));
    }

    #[test]
    fn pub_unused_report_parse_args_rejects_unknown_flag() {
        assert_eq!(None, parse_args(os_args(&["/ws/Cargo.toml", "--exclude", "a"])));
    }

    #[test]
    fn pub_unused_report_line_col_counts_characters() {
        let src = "fn a() {}\n// é\npub fn b() {}\n".as_bytes();
        let offset = "fn a() {}\n// é\npub fn ".len();
        assert_eq!(Some((3, 8)), line_col(src, offset));
    }

    #[test]
    fn pub_unused_report_line_col_rejects_offset_past_end() {
        assert_eq!(None, line_col(b"fn a() {}", 100));
    }

    fn record(
        src_path: &str,
        crate_name: &str,
        executable: bool,
        defs: &[Def],
        local_refs: &[&str],
        foreign_refs: &[&str],
    ) -> CrateRecord {
        CrateRecord {
            src_path: PathBuf::from(src_path),
            crate_name: crate_name.to_owned(),
            executable,
            defs: defs.to_vec(),
            local_refs: local_refs.iter().map(|key| (*key).to_owned()).collect(),
            foreign_refs: foreign_refs.iter().map(|key| (*key).to_owned()).collect(),
        }
    }

    fn def(key: &str, file: &str, lo: u32, parent: Option<&str>) -> Def {
        Def {
            key: key.to_owned(),
            file: PathBuf::from(file),
            lo,
            hi: lo + 1,
            descr: "function".to_owned(),
            path: key.to_owned(),
            parent: parent.map(str::to_owned),
        }
    }

    fn target(src_path: &str, required: bool) -> WorkspaceTarget {
        WorkspaceTarget {
            package: "p".to_owned(),
            name: "t".to_owned(),
            kind: "lib".to_owned(),
            src_path: PathBuf::from(src_path),
            required,
        }
    }

    fn keys<'a>(defs: &[&'a Def]) -> Vec<&'a str> {
        defs.iter().map(|def| def.key.as_str()).collect()
    }

    fn os_args(args: &[&str]) -> impl Iterator<Item = OsString> {
        args.iter().map(OsString::from).collect::<Vec<_>>().into_iter()
    }
}
