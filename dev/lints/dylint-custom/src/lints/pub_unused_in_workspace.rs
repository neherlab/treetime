//! Collects, per workspace crate, the reachable `pub` items it defines and the
//! workspace items it uses. The findings come from the `pub-unused-report`
//! binary, which runs after cargo has compiled every crate: a crate cannot
//! judge another crate's items while cargo may still be compiling their users.
//! Adapted from mordant's `unused_pub` (github.com/scarletindustries/mordant).

pub mod record;

use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

use rustc_hir::def::{DefKind, Res};
use rustc_hir::def_id::{DefId, LOCAL_CRATE, LocalDefId};
use rustc_hir::{
    Expr, ExprKind, HirId, ImplItem, ImplItemKind, Item, ItemKind, Node, Pat, PatExpr, PatExprKind,
    PatKind, Path as HirPath, QPath, TraitItem, TraitItemKind,
};
use rustc_lint::{LateContext, LateLintPass, LintContext as _};
use rustc_middle::middle::codegen_fn_attrs::CodegenFnAttrFlags;
use rustc_middle::ty::TyCtxt;
use rustc_middle::ty::print::with_no_trimmed_paths;
use rustc_session::config::CrateType;
use rustc_span::{FileName, Span};

use record::{CrateRecord, Def, RECORD_EXTENSION, RECORDS_DIR_ENV};

rustc_session::declare_lint! {
    /// Finds a reachable `pub` item that no crate of the workspace uses: no
    /// lib, bin, example, bench, or build script names it, calls it, or
    /// imports it, its own crate included. rustc's `dead_code` never reports
    /// such an item, because `pub` makes it reachable from outside the crate.
    /// Test builds take no part, so an item only tests use is reported.
    ///
    /// Not reported: items with `#[no_mangle]`, `#[export_name]` or `#[used]`,
    /// language items, the entry point, items produced by macros, items in a
    /// file brought in with `include!`, and items of proc-macro crates.
    pub PUB_UNUSED_IN_WORKSPACE,
    Warn,
    "a public item that no crate in the workspace uses"
}

pub struct PubUnusedInWorkspace {
    /// `None` when collection is off or the crate is a test build.
    dir: Option<PathBuf>,
    is_executable: bool,
    is_proc_macro: bool,
    defs: Vec<(LocalDefId, Def)>,
    local_refs: HashSet<LocalDefId>,
    foreign_refs: BTreeSet<String>,
    /// Crates with a record present when this crate started. cargo compiles
    /// every dependency first, so these include every workspace crate this
    /// crate can use.
    recorded_crates: HashSet<String>,
    keys: HashMap<DefId, String>,
}

impl PubUnusedInWorkspace {
    pub fn new() -> Self {
        Self {
            dir: None,
            is_executable: false,
            is_proc_macro: false,
            defs: Vec::new(),
            local_refs: HashSet::new(),
            foreign_refs: BTreeSet::new(),
            recorded_crates: HashSet::new(),
            keys: HashMap::new(),
        }
    }

    fn record_def(&mut self, cx: &LateContext<'_>, def_id: LocalDefId, item_span: Span, name: Span) {
        if self.dir.is_none() || self.is_proc_macro {
            return;
        }
        let hir_id = cx.tcx.local_def_id_to_hir_id(def_id);
        if item_span.from_expansion()
            || !cx.effective_visibilities.is_reachable(def_id)
            || exempt(cx, def_id)
            || included(cx, hir_id, item_span)
        {
            return;
        }
        let spec = cx.tcx.lint_level_spec_at_node(PUB_UNUSED_IN_WORKSPACE, hir_id);
        if spec.is_allow() || spec.is_expect() {
            if let Some(expectation) = spec.lint_id() {
                cx.fulfill_expectation(expectation);
            }
            return;
        }
        let Some((file, lo, hi)) = locate(cx, name) else {
            return;
        };
        let did = def_id.to_def_id();
        let crate_name = cx.tcx.crate_name(LOCAL_CRATE);
        let mut path = with_no_trimmed_paths!(cx.tcx.def_path_str(did));
        if !path.starts_with(&format!("{crate_name}::")) {
            path = format!("{crate_name}::{path}");
        }
        let def = Def {
            key: key(cx.tcx, did),
            file,
            lo,
            hi,
            descr: cx.tcx.def_descr(did).to_owned(),
            path,
            parent: parent_key(cx, did),
        };
        self.defs.push((def_id, def));
    }

    fn record_ref(&mut self, cx: &LateContext<'_>, def_id: DefId, from: HirId) {
        if self.dir.is_none() || from.owner.to_def_id() == def_id {
            return;
        }
        if let Some(local) = def_id.as_local() {
            self.local_refs.insert(local);
            return;
        }
        if !self
            .recorded_crates
            .contains(cx.tcx.crate_name(def_id.krate).as_str())
        {
            return;
        }
        let key = self.keys.entry(def_id).or_insert_with(|| key(cx.tcx, def_id));
        self.foreign_refs.insert(key.clone());
    }
}

rustc_session::impl_lint_pass!(PubUnusedInWorkspace => [PUB_UNUSED_IN_WORKSPACE]);

impl<'tcx> LateLintPass<'tcx> for PubUnusedInWorkspace {
    fn check_crate(&mut self, cx: &LateContext<'tcx>) {
        if cx.tcx.sess.is_test_crate() {
            return;
        }
        let Some(dir) = std::env::var_os(RECORDS_DIR_ENV).map(PathBuf::from) else {
            return;
        };
        let crate_types = cx.tcx.crate_types();
        self.is_executable = crate_types.contains(&CrateType::Executable);
        self.is_proc_macro = crate_types.contains(&CrateType::ProcMacro);
        self.recorded_crates = recorded_crates(&dir);
        self.dir = Some(dir);
    }

    fn check_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx Item<'tcx>) {
        let checked = matches!(
            item.kind,
            ItemKind::Fn { .. }
                | ItemKind::Struct(..)
                | ItemKind::Enum(..)
                | ItemKind::Union(..)
                | ItemKind::Const(..)
                | ItemKind::Static(..)
                | ItemKind::TyAlias(..)
                | ItemKind::Trait { .. }
        );
        if checked && let Some(ident) = item.kind.ident() {
            self.record_def(cx, item.owner_id.def_id, item.span, ident.span);
        }
    }

    fn check_impl_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx ImplItem<'tcx>) {
        if matches!(item.kind, ImplItemKind::Fn(..) | ImplItemKind::Const(..))
            && let Node::Item(parent) = cx.tcx.parent_hir_node(item.hir_id())
            && let ItemKind::Impl(imp) = parent.kind
            && imp.of_trait.is_none()
        {
            self.record_def(cx, item.owner_id.def_id, item.span, item.ident.span);
        }
    }

    fn check_trait_item(&mut self, cx: &LateContext<'tcx>, item: &'tcx TraitItem<'tcx>) {
        if matches!(item.kind, TraitItemKind::Fn(..) | TraitItemKind::Const(..)) {
            self.record_def(cx, item.owner_id.def_id, item.span, item.ident.span);
        }
    }

    fn check_path(&mut self, cx: &LateContext<'tcx>, path: &HirPath<'tcx>, hir_id: HirId) {
        let Res::Def(_, def_id) = path.res else {
            return;
        };
        // `impl Foo { .. }` and `impl Trait for Foo` do not use `Foo`. They do
        // use an alias written there: the impl is of the type it names.
        if !matches!(path.res, Res::Def(DefKind::TyAlias, _))
            && let Node::Item(item) = cx.tcx.parent_hir_node(hir_id)
            && let ItemKind::Impl(imp) = item.kind
            && imp.self_ty.hir_id == hir_id
        {
            return;
        }
        let def_id = owning_item(cx, def_id);
        if path.span.in_derive_expansion() && enclosing_impl_self_adt(cx, hir_id) == Some(def_id) {
            return;
        }
        self.record_ref(cx, def_id, hir_id);
    }

    fn check_expr(&mut self, cx: &LateContext<'tcx>, expr: &'tcx Expr<'tcx>) {
        let def_id = match expr.kind {
            ExprKind::MethodCall(..) => cx.typeck_results().type_dependent_def_id(expr.hir_id),
            ExprKind::Path(ref qpath @ QPath::TypeRelative(..))
            | ExprKind::Struct(&ref qpath @ QPath::TypeRelative(..), ..) => {
                cx.qpath_res(qpath, expr.hir_id).opt_def_id()
            },
            _ => None,
        };
        if let Some(def_id) = def_id {
            self.record_ref(cx, owning_item(cx, def_id), expr.hir_id);
        }
    }

    fn check_pat(&mut self, cx: &LateContext<'tcx>, pat: &'tcx Pat<'tcx>) {
        let (qpath, hir_id) = match pat.kind {
            PatKind::Struct(ref qpath @ QPath::TypeRelative(..), ..)
            | PatKind::TupleStruct(ref qpath @ QPath::TypeRelative(..), ..) => (qpath, pat.hir_id),
            PatKind::Expr(PatExpr {
                hir_id,
                kind: PatExprKind::Path(qpath @ QPath::TypeRelative(..)),
                ..
            }) => (qpath, *hir_id),
            _ => return,
        };
        if let Some(def_id) = cx.qpath_res(qpath, hir_id).opt_def_id() {
            self.record_ref(cx, owning_item(cx, def_id), hir_id);
        }
    }

    fn check_crate_post(&mut self, cx: &LateContext<'tcx>) {
        let Some(dir) = self.dir.take() else {
            return;
        };
        let Some(src_path) = crate_root(cx) else {
            return;
        };
        let defs = std::mem::take(&mut self.defs);
        let local_refs = defs
            .iter()
            .filter(|(id, _)| self.local_refs.contains(id))
            .map(|(_, def)| def.key.clone())
            .collect();
        let record = CrateRecord {
            src_path,
            crate_name: cx.tcx.crate_name(LOCAL_CRATE).to_string(),
            executable: self.is_executable,
            defs: defs.into_iter().map(|(_, def)| def).collect(),
            local_refs,
            foreign_refs: std::mem::take(&mut self.foreign_refs),
        };
        if let Err(err) = write_record(&dir, &record) {
            cx.sess().dcx().err(format!(
                "pub_unused_in_workspace: cannot write the record of `{}` under `{}`: {err}",
                record.crate_name,
                dir.display()
            ));
        }
    }
}

/// Crate names of the library records in `dir`. A record that cannot be read
/// is skipped here; the report binary names it.
fn recorded_crates(dir: &Path) -> HashSet<String> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return HashSet::new();
    };
    entries
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|ext| ext == RECORD_EXTENSION))
        .filter_map(|path| std::fs::read(path).ok())
        .filter_map(|bytes| serde_json::from_slice::<CrateRecord>(&bytes).ok())
        .filter(|record| !record.executable)
        .map(|record| record.crate_name)
        .collect()
}

/// Written under a temporary name and renamed, so a reader never sees half a
/// record.
fn write_record(dir: &Path, record: &CrateRecord) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)?;
    let path = record_path(dir, &record.src_path);
    let tmp = path.with_extension(format!("tmp{}", std::process::id()));
    std::fs::write(&tmp, serde_json::to_vec(record).map_err(std::io::Error::other)?)?;
    std::fs::rename(&tmp, &path)
}

/// Record file of the crate rooted at `src_path`: the path escaped byte by
/// byte, so distinct roots never share a file.
pub fn record_path(dir: &Path, src_path: &Path) -> PathBuf {
    let tag = src_path
        .as_os_str()
        .as_encoded_bytes()
        .iter()
        .map(|&byte| {
            if byte.is_ascii_alphanumeric() {
                char::from(byte).to_string()
            } else {
                format!("_{byte:02x}")
            }
        })
        .collect::<String>();
    dir.join(format!("{tag}.{RECORD_EXTENSION}"))
}

/// Absolute path of the crate root file, as `cargo metadata` reports it.
fn crate_root(cx: &LateContext<'_>) -> Option<PathBuf> {
    let root = cx.tcx.sess.local_crate_source_file()?;
    absolute(root.local_path()?)
}

fn absolute(path: &Path) -> Option<PathBuf> {
    if path.is_absolute() {
        return Some(path.to_path_buf());
    }
    Some(std::env::current_dir().ok()?.join(path))
}

/// Absolute file and byte range of `span`, if it is in a real file.
fn locate(cx: &LateContext<'_>, span: Span) -> Option<(PathBuf, u32, u32)> {
    let sm = cx.tcx.sess.source_map();
    let file = sm.lookup_source_file(span.lo());
    let FileName::Real(real) = &file.name else {
        return None;
    };
    Some((
        absolute(real.local_path()?)?,
        (span.lo() - file.start_pos).0,
        (span.hi() - file.start_pos).0,
    ))
}

/// Crate name plus definition path: the same string whichever crate computes it.
fn key(tcx: TyCtxt<'_>, def_id: DefId) -> String {
    format!(
        "{}{}",
        tcx.crate_name(def_id.krate),
        tcx.def_path(def_id).to_string_no_crate_verbose()
    )
}

/// The item a use counts for: a constructor or variant counts as its struct or
/// enum.
fn owning_item(cx: &LateContext<'_>, def_id: DefId) -> DefId {
    match cx.tcx.def_kind(def_id) {
        DefKind::Ctor(..) => owning_item(cx, cx.tcx.parent(def_id)),
        DefKind::Variant => cx.tcx.parent(def_id),
        _ => def_id,
    }
}

/// The trait or, for an inherent impl, the type an associated item belongs to.
fn parent_key(cx: &LateContext<'_>, did: DefId) -> Option<String> {
    if !matches!(cx.tcx.def_kind(did), DefKind::AssocFn | DefKind::AssocConst { .. }) {
        return None;
    }
    let parent = cx.tcx.parent(did);
    let owner = match cx.tcx.def_kind(parent) {
        DefKind::Trait => parent,
        DefKind::Impl { .. } => cx
            .tcx
            .type_of(parent)
            .instantiate_identity()
            .skip_normalization()
            .ty_adt_def()?
            .did(),
        _ => return None,
    };
    Some(key(cx.tcx, owner))
}

/// The ADT of the `impl` block `hir_id` is written in, if any.
fn enclosing_impl_self_adt(cx: &LateContext<'_>, hir_id: HirId) -> Option<DefId> {
    let owner = hir_id.owner.to_def_id();
    let impl_id = match cx.tcx.def_kind(owner) {
        DefKind::Impl { .. } => owner,
        DefKind::AssocFn | DefKind::AssocConst { .. } | DefKind::AssocTy => cx.tcx.parent(owner),
        _ => return None,
    };
    if !matches!(cx.tcx.def_kind(impl_id), DefKind::Impl { .. }) {
        return None;
    }
    cx.tcx
        .type_of(impl_id)
        .instantiate_identity()
        .skip_normalization()
        .ty_adt_def()
        .map(|adt| adt.did())
}

/// In a different file from the module body around it: brought in by
/// `include!`, which is how generated code enters a crate. An out-of-line
/// `mod m;` also changes file, but then the module body and its items agree.
fn included(cx: &LateContext<'_>, hir_id: HirId, item_span: Span) -> bool {
    let sm = cx.tcx.sess.source_map();
    let file_of = |span: Span| sm.lookup_source_file(span.lo()).start_pos;
    let mut file = file_of(item_span);
    let mut module = cx.tcx.parent_module(hir_id);
    loop {
        match cx.tcx.hir_node_by_def_id(module.to_local_def_id()) {
            Node::Crate(body) => return file != file_of(body.spans.inner_span),
            Node::Item(item) if let ItemKind::Mod(_, body) = item.kind => {
                if file != file_of(body.spans.inner_span) {
                    return true;
                }
                file = file_of(item.span);
                module = cx.tcx.parent_module(item.hir_id());
            },
            _ => return false,
        }
    }
}

/// Reached by symbol name from outside Rust, required by the language, or the
/// program's entry point.
fn exempt(cx: &LateContext<'_>, def_id: LocalDefId) -> bool {
    let did = def_id.to_def_id();
    let by_symbol = matches!(
        cx.tcx.def_kind(did),
        DefKind::Fn | DefKind::AssocFn | DefKind::Static { .. }
    ) && {
        let attrs = cx.tcx.codegen_fn_attrs(did);
        attrs.contains_extern_indicator()
            || attrs
                .flags
                .intersects(CodegenFnAttrFlags::USED_COMPILER | CodegenFnAttrFlags::USED_LINKER)
    };
    by_symbol
        || cx.tcx.lang_items().from_def_id(did).is_some()
        || cx.tcx.entry_fn(()).is_some_and(|(entry, _)| entry == did)
}
