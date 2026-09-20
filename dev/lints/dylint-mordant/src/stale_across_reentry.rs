use clippy_utils::res::MaybeResPath;
use clippy_utils::visitors::{for_each_expr, for_each_expr_without_closures, is_local_used};
use rustc_hir::def_id::DefId;
use rustc_hir::intravisit::{Visitor, walk_expr};
use rustc_hir::{
    BinOpKind, Block, Expr, ExprKind, HirId, ImplicitSelfKind, MatchSource, Mutability, PatKind,
    Stmt, StmtKind, UnOp,
};
use rustc_lint::{LateContext, LateLintPass};
use rustc_middle::ty;
use rustc_span::Span;
use rustc_span::symbol::{Symbol, sym};

use crate::MordantConfig;
use crate::adt_facts::impl_self_adt;
use crate::baseline::emit_with_note;
use crate::hir_shapes::{
    FieldChain, def_path_names, dotted, field_chain, is_self_path, stmt_expr,
    strip_generic_segments,
};

rustc_session::declare_lint! {
    /// Flags a value read from a field of `self` (`let n = self.items.len()`,
    /// `let p = self.buf.as_ptr()`, `let had = self.cb.is_some()`), then a
    /// call that can run other code with access to `self`, then a later
    /// statement that uses the value as if the field had not changed:
    /// indexing or removing at `n`, unwrapping under `had`, or touching `p`
    /// at all. Whatever ran in between may have pushed, popped, taken, or
    /// reallocated. The panic or dangling access only shows up on the
    /// re-entrant path, which is the one tests skip.
    ///
    /// The re-entrant statements are calls through a closure or fn-pointer
    /// parameter or field, calls on a `dyn Trait` receiver, `.await`, and
    /// anything named in `stale-across-reentry-callees` (a project's own
    /// dispatch points, such as the methods that run script callbacks). A
    /// plain call to a named function is not one: whether it re-enters is not
    /// decidable from here, so it is only counted when the config says so.
    /// A configured callee counts wherever it appears; the language-level
    /// ones count only where the language lets them reach the field: in a
    /// `&mut self` method the call must actually receive `self` or the field
    /// by `&mut` (an argument coerced to `&Host` or `&Vec<_>` does not
    /// count), and in a `&self` method, or through a shared reference, the
    /// field's type must be changeable behind `&`: a type of this crate that
    /// is not `Freeze` or reaches non-`Freeze` state through a pointer. A
    /// std collection, or a `Freeze` newtype over one, changes only through
    /// `&mut`. A call counts only on a path that goes on to the later
    /// statement: one inside a branch that returns, breaks out of the block,
    /// or panics does not, and one inside a closure that is merely built
    /// here has not run.
    ///
    /// Silent when the fact is about a local rather than a field of `self`,
    /// when the field is a reference or a pointer wrapper rather than an
    /// owned collection (`NonNull::as_ptr` copies an address out, it does
    /// not point into the field), when the field is re-queried after the
    /// re-entry (a fresh `len()`, `is_some()`, or a whole-field assignment)
    /// before the fact is reused, when the reuse cannot panic (`get(n)`,
    /// `truncate(n)`), and when the fact and its reuse are not in the same
    /// block.
    pub STALE_ACROSS_REENTRY,
    Warn,
    "a fact about a field of self reused after a call that can re-enter and change it"
}

pub struct StaleAcrossReentry {
    pub config: &'static MordantConfig,
}

rustc_session::impl_lint_pass!(StaleAcrossReentry => [STALE_ACROSS_REENTRY]);

impl StaleAcrossReentry {
    /// Whether the call to `def` with `args` is one the config names, under
    /// any spelling of either the item typeck resolved or the impl item the
    /// receiver type sends it to.
    fn configured<'tcx>(
        &self,
        cx: &LateContext<'tcx>,
        def: DefId,
        args: ty::GenericArgsRef<'tcx>,
    ) -> bool {
        let callees = &self.config.stale_across_reentry_callees;
        if callees.is_empty() {
            return false;
        }
        std::iter::once(def)
            .chain(impl_item_of(cx, def, args))
            .flat_map(|def| callee_names(cx, def))
            .any(|[name, with_crate]| {
                callees
                    .iter()
                    .any(|p| matches_pattern(&name, p) || matches_pattern(&with_crate, p))
            })
    }
}

/// A call to a trait method is recorded against the trait's item, whose path
/// is `Runner::run_job` whatever the receiver; when the receiver type picks
/// the impl here, this is the impl's item, so that `Worker::run_job` can be
/// named too. Dispatch through `dyn` or a type parameter has no impl to name.
fn impl_item_of<'tcx>(
    cx: &LateContext<'tcx>,
    def: DefId,
    args: ty::GenericArgsRef<'tcx>,
) -> Option<DefId> {
    cx.tcx.trait_of_assoc(def)?;
    let args = cx
        .tcx
        .try_normalize_erasing_regions(cx.typing_env(), ty::Unnormalized::new_wip(args))
        .ok()?;
    let instance = ty::Instance::try_resolve(cx.tcx, cx.typing_env(), def, args).ok()??;
    match instance.def {
        ty::InstanceKind::Item(item) if item != def => Some(item),
        _ => None,
    }
}

/// The spellings a pattern may match a callee under, each also prefixed
/// with its crate name. The def path covers free functions, inherent methods
/// (`Vm::run_callback`) and trait items (`Runner::run_job`); an impl's item
/// renders as `<Worker as Runner>::run_job`, which no `Type::method` pattern
/// matches, so it is spelled `Worker::run_job` and `Runner::run_job` instead.
fn callee_names(cx: &LateContext<'_>, def: DefId) -> Vec<[String; 2]> {
    let mut names = vec![def_path_names(cx, def)];
    if let Some(impl_did) = cx.tcx.impl_of_assoc(def)
        && let Some(trait_ref) = cx.tcx.impl_opt_trait_ref(impl_did)
    {
        let (krate, item) = (cx.tcx.crate_name(def.krate), cx.tcx.item_name(def));
        let owners = impl_self_adt(cx, impl_did)
            .map(|adt| adt.did())
            .into_iter()
            .chain([trait_ref.skip_binder().def_id]);
        names.extend(owners.map(|owner| {
            let name = format!(
                "{}::{item}",
                strip_generic_segments(&cx.tcx.def_path_str(owner))
            );
            let with_crate = format!("{krate}::{name}");
            [name, with_crate]
        }));
    }
    names
}

/// Segment-wise suffix match: `run_callback` and `Vm::run_callback` both
/// match `bun::vm::Vm::run_callback`; a final `dispatch*` matches any last
/// segment starting with `dispatch`. Substrings never match by accident.
fn matches_pattern(name: &str, pattern: &str) -> bool {
    let mut want = pattern.rsplit("::");
    let mut have = name.rsplit("::");
    let (Some(want_last), Some(have_last)) = (want.next(), have.next()) else {
        return false;
    };
    let last_matches = match want_last.strip_suffix('*') {
        Some(prefix) => have_last.starts_with(prefix),
        None => have_last == want_last,
    };
    last_matches && want.all(|w| have.next() == Some(w))
}

/// What kind of fact a binding holds about its place, which decides what
/// counts as reusing it.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Fact {
    /// `len()` / `capacity()`: reused by indexing or a position-taking
    /// method on the place, or by gating one of those on the number.
    Count,
    /// `is_empty()` / `is_some()` / `is_none()`: reused by gating an access
    /// to the place on the flag.
    Flag,
    /// `as_ptr()` / `as_mut_ptr()`: reused by any mention at all, since the
    /// allocation it points into is what a re-entrant callee can move.
    Pointer,
}

impl Fact {
    fn of(method: &str) -> Option<Self> {
        match method {
            "len" | "capacity" => Some(Fact::Count),
            "is_empty" | "is_some" | "is_none" => Some(Fact::Flag),
            "as_ptr" | "as_mut_ptr" => Some(Fact::Pointer),
            _ => None,
        }
    }

    /// The binding's type must be what the method name promises; a project
    /// method that happens to share the name but returns something richer is
    /// not one of these facts.
    fn matches_ty(self, ty: ty::Ty<'_>) -> bool {
        match self {
            Fact::Count => ty.is_integral(),
            Fact::Flag => ty.is_bool(),
            Fact::Pointer => ty.is_raw_ptr(),
        }
    }

    fn help(self, name: Symbol, place: &str) -> String {
        match self {
            Fact::Count | Fact::Flag => format!(
                "read `{place}` again after the call instead of reusing `{name}`, or hold the element itself across the call"
            ),
            Fact::Pointer => format!(
                "take the pointer from `{place}` again after the call, or keep the data alive across it by value"
            ),
        }
    }
}

/// Methods that take a position or otherwise panic / are unsafe when the
/// container has shrunk since the position was computed. Non-panicking
/// lookups (`get`, `truncate`) are deliberately absent.
const POSITIONAL: &[&str] = &[
    "remove",
    "swap_remove",
    "insert",
    "split_off",
    "drain",
    "split_at",
    "split_at_mut",
    "swap",
    "copy_within",
    "set_len",
    "get_unchecked",
    "get_unchecked_mut",
];

const UNWRAPS: &[&str] = &["unwrap", "expect", "unwrap_unchecked"];

/// Adapters between an `Option`/`Vec` field and the `unwrap` that bets on
/// it: `self.cb.as_ref().unwrap()`, `self.items.get(n).unwrap()`.
const ADAPTERS: &[&str] = &[
    "as_ref",
    "as_mut",
    "as_deref",
    "as_deref_mut",
    "take",
    "get",
    "get_mut",
    "first",
    "last",
    "first_mut",
    "last_mut",
];

/// `self.a.b` as `[a, b]`; anything not rooted at `self` through fields
/// alone has no identity a re-entrant callee shares with this function.
fn self_place(e: &Expr<'_>) -> Option<Vec<Symbol>> {
    let FieldChain { root, fields } = field_chain(e);
    (is_self_path(root) && !fields.is_empty()).then_some(fields)
}

/// What the fact's field is, which decides who could change it underneath
/// the fact.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Container {
    /// A std collection or `Option`, or a type of this crate with nothing
    /// mutable behind `&`. Changes only through `&mut`, so a callee has to be
    /// handed one.
    Frozen,
    /// A type of this crate that can change behind `&`, so anyone holding a
    /// shared reference to the owner can change it.
    Interior,
}

/// The field type a fact can be about. A reference-typed field (`&[u8]`)
/// cannot change under a borrow at all; a pointer wrapper's `as_ptr()`
/// (`NonNull`, `RefPtr`) copies a stored address out rather than pointing
/// into storage the field owns, so a pointer fact needs a `Vec` or `String`.
fn container_of<'tcx>(cx: &LateContext<'tcx>, ty: ty::Ty<'tcx>, fact: Fact) -> Option<Container> {
    let ty::Adt(adt, _) = ty.kind() else {
        return None;
    };
    if adt.did().is_local() {
        if fact == Fact::Pointer {
            return None;
        }
        return Some(if changes_behind_shared(cx, ty, 4) {
            Container::Interior
        } else {
            Container::Frozen
        });
    }
    let name = cx.tcx.get_diagnostic_name(adt.did())?;
    let owning = match fact {
        Fact::Pointer => matches!(name.as_str(), "Vec" | "String"),
        Fact::Count | Fact::Flag => matches!(
            name.as_str(),
            "Vec"
                | "String"
                | "VecDeque"
                | "HashMap"
                | "HashSet"
                | "BTreeMap"
                | "BTreeSet"
                | "Option"
        ),
    };
    owning.then_some(Container::Frozen)
}

/// Whether a `&` to the type reaches something that can change: the type
/// itself is not `Freeze`, or one of its fields gets to non-`Freeze` state
/// through a reference, a raw pointer, or a pointer type's argument
/// (`Box<RefCell<_>>`, `Rc<Cell<_>>`), which `Freeze` does not look
/// through. A bare type parameter met behind a pointer is unknown rather
/// than mutable; `depth` bounds the walk through self-referential types.
fn changes_behind_shared<'tcx>(cx: &LateContext<'tcx>, ty: ty::Ty<'tcx>, depth: usize) -> bool {
    if depth == 0 {
        return false;
    }
    let inner = |t| changes_behind_shared(cx, t, depth - 1);
    match ty.kind() {
        ty::Param(_) => false,
        _ if !ty.is_freeze(cx.tcx, cx.typing_env()) => true,
        ty::RawPtr(..) => true,
        ty::Ref(_, pointee, _) | ty::Array(pointee, _) | ty::Slice(pointee) => inner(*pointee),
        ty::Tuple(elems) => elems.iter().any(inner),
        ty::Adt(adt, args) if adt.did().is_local() => adt
            .all_fields()
            .any(|f| inner(f.ty(cx.tcx, args).skip_normalization())),
        ty::Adt(adt, _) if adt.is_phantom_data() => false,
        ty::Adt(adt, args) => {
            cx.tcx.is_diagnostic_item(sym::NonNull, adt.did()) || args.types().any(inner)
        }
        _ => false,
    }
}

struct Tracked {
    binding: HirId,
    name: Symbol,
    place: Vec<Symbol>,
    fact: Fact,
    container: Container,
}

/// The fact query inside a `let` initializer, looking through the wrappers
/// facts are habitually stored under: `as u32`, `!`, `- 1`.
fn fact_query<'tcx>(
    cx: &LateContext<'tcx>,
    e: &'tcx Expr<'tcx>,
) -> Option<(Vec<Symbol>, Fact, Container)> {
    match &e.kind {
        ExprKind::MethodCall(seg, recv, [], _) => {
            let fact = Fact::of(seg.ident.as_str())?;
            if !fact.matches_ty(cx.typeck_results().expr_ty(e)) {
                return None;
            }
            let container = container_of(cx, cx.typeck_results().expr_ty(recv), fact)?;
            Some((self_place(recv)?, fact, container))
        }
        ExprKind::Cast(inner, _) | ExprKind::Unary(UnOp::Not, inner) => fact_query(cx, inner),
        ExprKind::Binary(op, l, r)
            if matches!(
                op.node,
                BinOpKind::Add | BinOpKind::Sub | BinOpKind::Mul | BinOpKind::Div
            ) =>
        {
            fact_query(cx, l).or_else(|| fact_query(cx, r))
        }
        ExprKind::DropTemps(inner) => fact_query(cx, inner),
        _ => None,
    }
}

fn tracked_of<'tcx>(cx: &LateContext<'tcx>, stmt: &'tcx Stmt<'tcx>) -> Option<Tracked> {
    let StmtKind::Let(l) = stmt.kind else {
        return None;
    };
    let PatKind::Binding(_, binding, ident, None) = l.pat.kind else {
        return None;
    };
    let (place, fact, container) = fact_query(cx, l.init?)?;
    // A derived value must still be the fact's own kind: `len() > 0` is a
    // Flag spelled through Count, and gating on it is handled by tracking
    // nothing rather than by guessing.
    if !fact.matches_ty(cx.typeck_results().pat_ty(l.pat)) {
        return None;
    }
    Some(Tracked {
        binding,
        name: ident.name,
        place,
        fact,
        container,
    })
}

/// How the enclosing method holds `self`, which is what bounds who else can
/// get at its fields while the method runs.
#[derive(Clone, Copy, PartialEq, Eq)]
enum SelfAccess {
    /// `&mut self` or `self`: nothing else can reach the fields unless this
    /// method hands them out.
    Exclusive,
    /// `&self`: every other holder of a reference can reach them too, and
    /// changes them if the field's type allows it.
    Shared,
}

fn self_access(cx: &LateContext<'_>, block: &Block<'_>) -> Option<SelfAccess> {
    let owner = cx.tcx.hir_get_parent_item(block.hir_id);
    match cx.tcx.hir_owner_node(owner).fn_decl()?.implicit_self() {
        ImplicitSelfKind::Imm | ImplicitSelfKind::Mut | ImplicitSelfKind::RefMut => {
            Some(SelfAccess::Exclusive)
        }
        ImplicitSelfKind::RefImm => Some(SelfAccess::Shared),
        ImplicitSelfKind::None => None,
    }
}

/// `self`, `&mut self.a`, `&mut self.a.b` handed to a call, where the
/// tracked place is `self.a.b`: the callee can now change it. What the
/// callee receives is read off the argument's adjusted type, so `self`
/// passed to an `fn(&Host)`, or `&mut self.a` to an `fn(&A)`, arrives shared
/// and counts only when the field's type can be changed through a shared
/// reference.
fn hands_out<'tcx>(cx: &LateContext<'tcx>, arg: &'tcx Expr<'tcx>, t: &Tracked) -> bool {
    let mut e = arg;
    let mut shared = matches!(
        cx.typeck_results().expr_ty_adjusted(arg).kind(),
        ty::Ref(_, _, Mutability::Not)
    );
    let counts = |shared: bool| !shared || t.container == Container::Interior;
    loop {
        match &e.kind {
            ExprKind::AddrOf(_, Mutability::Mut, inner) => {
                if self_place(inner).is_some_and(|p| t.place.starts_with(&p)) {
                    return counts(shared);
                }
                e = inner;
            }
            ExprKind::AddrOf(_, Mutability::Not, inner) => {
                shared = true;
                e = inner;
            }
            ExprKind::Unary(UnOp::Deref, inner) | ExprKind::DropTemps(inner) => e = inner,
            _ => return is_self_path(e) && counts(shared),
        }
    }
}

/// How an expression gives control away.
enum Exit<'tcx> {
    /// A configured callee: the project says it re-enters, and re-entry
    /// reaches whatever the project's other handles reach, so no argument
    /// analysis applies.
    Configured,
    /// Code this crate cannot see (a closure or fn pointer, a `dyn` method,
    /// the executor behind `.await`), which can only touch what the language
    /// lets it: what it is handed, plus whatever a shared reference reaches.
    Opaque { args: &'tcx [Expr<'tcx>] },
}

impl StaleAcrossReentry {
    fn exit_of<'tcx>(&self, cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>) -> Option<Exit<'tcx>> {
        match &e.kind {
            ExprKind::Call(callee, args) => {
                let ty = cx.typeck_results().expr_ty_adjusted(callee).peel_refs();
                let ty = ty.boxed_ty().unwrap_or(ty).peel_refs();
                match ty.kind() {
                    ty::FnPtr(..) | ty::Dynamic(..) | ty::Param(_) => Some(Exit::Opaque { args }),
                    ty::FnDef(def, fn_args) if self.configured(cx, *def, fn_args) => {
                        Some(Exit::Configured)
                    }
                    _ => None,
                }
            }
            ExprKind::MethodCall(_, recv, args, _) => {
                if cx
                    .typeck_results()
                    .type_dependent_def_id(e.hir_id)
                    .is_some_and(|def| {
                        self.configured(cx, def, cx.typeck_results().node_args(e.hir_id))
                    })
                {
                    Some(Exit::Configured)
                } else if cx
                    .typeck_results()
                    .expr_ty_adjusted(recv)
                    .peel_refs()
                    .is_trait()
                {
                    Some(Exit::Opaque { args })
                } else {
                    None
                }
            }
            ExprKind::Match(_, _, MatchSource::AwaitDesugar) => Some(Exit::Opaque { args: &[] }),
            _ => None,
        }
    }

    fn reaches<'tcx>(
        &self,
        cx: &LateContext<'tcx>,
        e: &'tcx Expr<'tcx>,
        t: &Tracked,
        access: SelfAccess,
    ) -> bool {
        match self.exit_of(cx, e) {
            None => false,
            Some(Exit::Configured) => true,
            Some(Exit::Opaque { args }) => match access {
                SelfAccess::Exclusive => args.iter().any(|a| hands_out(cx, a, t)),
                SelfAccess::Shared => t.container == Container::Interior,
            },
        }
    }

    /// The first expression in the statement `e` that gives control to code
    /// which can change the tracked field, on a path that then goes on to
    /// the statements after `e`.
    fn reentry<'tcx>(
        &self,
        cx: &LateContext<'tcx>,
        e: &'tcx Expr<'tcx>,
        t: &Tracked,
        access: SelfAccess,
    ) -> Option<Span> {
        let mut v = Reentry {
            lint: self,
            cx,
            t,
            access,
            targets: Vec::new(),
            found: None,
        };
        v.visit_expr(e);
        v.found
    }
}

/// Walks one statement for a re-entry that falls through to the next. As a
/// HIR visitor it stays out of closure bodies, whose calls have not run. A
/// subtree of type `!` never completes; the statements after it are reached
/// from inside it only by a `break` or `continue` to a loop or labelled
/// block that is itself inside the statement (`targets`), so a diverging
/// subtree without one, a returning branch or a `?` failure arm, is skipped
/// along with any call inside it.
struct Reentry<'a, 'tcx> {
    lint: &'a StaleAcrossReentry,
    cx: &'a LateContext<'tcx>,
    t: &'a Tracked,
    access: SelfAccess,
    targets: Vec<HirId>,
    found: Option<Span>,
}

impl<'tcx> Reentry<'_, 'tcx> {
    fn diverges(&self, e: &'tcx Expr<'tcx>) -> bool {
        self.cx
            .typeck_results()
            .expr_ty_opt(e)
            .is_some_and(ty::Ty::is_never)
    }

    fn escapes_into_statement(&self, e: &'tcx Expr<'tcx>) -> bool {
        for_each_expr_without_closures(e, |inner: &Expr<'tcx>| {
            let target = match &inner.kind {
                ExprKind::Break(dest, _) | ExprKind::Continue(dest) => dest.target_id.ok(),
                _ => None,
            };
            if target.is_some_and(|id| self.targets.contains(&id)) {
                std::ops::ControlFlow::Break(())
            } else {
                std::ops::ControlFlow::Continue(())
            }
        })
        .is_some()
    }
}

impl<'tcx> Visitor<'tcx> for Reentry<'_, 'tcx> {
    fn visit_expr(&mut self, e: &'tcx Expr<'tcx>) {
        if self.found.is_some() || (self.diverges(e) && !self.escapes_into_statement(e)) {
            return;
        }
        if self.lint.reaches(self.cx, e, self.t, self.access) {
            self.found = Some(e.span);
            return;
        }
        match &e.kind {
            ExprKind::Loop(..) | ExprKind::Block(_, Some(_)) => {
                self.targets.push(e.hir_id);
                walk_expr(self, e);
                self.targets.pop();
            }
            _ => walk_expr(self, e),
        }
    }
}

/// The statement invalidates the tracked fact on its own account, before any
/// re-entry is involved: it assigns the binding, the place, or a prefix of
/// the place. After a re-entry a fresh query of the place (`after_reentry`)
/// counts too — the code re-derived the fact, and what follows is its own.
fn refreshes<'tcx>(
    cx: &LateContext<'tcx>,
    e: &'tcx Expr<'tcx>,
    t: &Tracked,
    after_reentry: bool,
) -> bool {
    for_each_expr(cx, e, |inner: &Expr<'tcx>| {
        let yes = match &inner.kind {
            ExprKind::MethodCall(seg, recv, [], _)
                if after_reentry && Fact::of(seg.ident.as_str()).is_some() =>
            {
                self_place(recv).is_some_and(|p| p == t.place)
            }
            ExprKind::Assign(lhs, _, _) | ExprKind::AssignOp(_, lhs, _) => {
                lhs.res_local_id() == Some(t.binding)
                    || self_place(lhs).is_some_and(|p| t.place.starts_with(&p))
            }
            _ => false,
        };
        if yes {
            std::ops::ControlFlow::Break(())
        } else {
            std::ops::ControlFlow::Continue(())
        }
    })
    .is_some()
}

struct Reuse<'a, 'tcx> {
    cx: &'a LateContext<'tcx>,
    t: &'a Tracked,
    /// Depth of enclosing branches whose condition read the binding.
    gated: usize,
    found: Option<Span>,
}

impl<'a, 'tcx> Reuse<'a, 'tcx> {
    fn depends(&self, args: &'tcx [Expr<'tcx>]) -> bool {
        self.gated > 0
            || args
                .iter()
                .any(|a| is_local_used(self.cx, a, self.t.binding))
    }

    fn access(&self, e: &'tcx Expr<'tcx>) -> bool {
        match &e.kind {
            ExprKind::Index(base, idx, _) => {
                self_place(base).is_some_and(|p| p == self.t.place)
                    && self.depends(std::slice::from_ref(*idx))
            }
            ExprKind::MethodCall(seg, recv, args, _) => {
                let name = seg.ident.as_str();
                if POSITIONAL.contains(&name) {
                    return self_place(recv).is_some_and(|p| p == self.t.place)
                        && self.depends(args);
                }
                if !UNWRAPS.contains(&name) {
                    return false;
                }
                let mut depends = self.gated > 0;
                let mut inner: &Expr<'tcx> = recv;
                while let ExprKind::MethodCall(s, r, a, _) = &inner.kind
                    && ADAPTERS.contains(&s.ident.as_str())
                {
                    depends |= self.depends(a);
                    inner = r;
                }
                depends && self_place(inner).is_some_and(|p| p == self.t.place)
            }
            _ => false,
        }
    }
}

impl<'tcx> Visitor<'tcx> for Reuse<'_, 'tcx> {
    fn visit_expr(&mut self, e: &'tcx Expr<'tcx>) {
        if self.found.is_some() {
            return;
        }
        if self.t.fact == Fact::Pointer {
            if e.res_local_id() == Some(self.t.binding) {
                self.found = Some(e.span);
            } else {
                walk_expr(self, e);
            }
            return;
        }
        if self.access(e) {
            self.found = Some(e.span);
            return;
        }
        // A condition that reads the number may access the field itself
        // (`if self.items[n] > 0`, `match self.items.remove(n)`); if it does
        // not, the branches under it are gated on the number.
        match &e.kind {
            ExprKind::If(cond, then, els) if is_local_used(self.cx, *cond, self.t.binding) => {
                self.visit_expr(cond);
                self.gated += 1;
                self.visit_expr(then);
                if let Some(els) = els {
                    self.visit_expr(els);
                }
                self.gated -= 1;
            }
            ExprKind::Match(scrut, arms, _) if is_local_used(self.cx, *scrut, self.t.binding) => {
                self.visit_expr(scrut);
                self.gated += 1;
                for arm in *arms {
                    self.visit_expr(arm.body);
                }
                self.gated -= 1;
            }
            _ => walk_expr(self, e),
        }
    }
}

fn reuse_in<'tcx>(cx: &LateContext<'tcx>, e: &'tcx Expr<'tcx>, t: &Tracked) -> Option<Span> {
    let mut v = Reuse {
        cx,
        t,
        gated: 0,
        found: None,
    };
    v.visit_expr(e);
    v.found
}

impl<'tcx> LateLintPass<'tcx> for StaleAcrossReentry {
    fn check_block(&mut self, cx: &LateContext<'tcx>, block: &'tcx Block<'tcx>) {
        let mut access = None;
        for (i, stmt) in block.stmts.iter().enumerate() {
            let Some(t) = tracked_of(cx, stmt) else {
                continue;
            };
            let Some(access) = *access.get_or_insert_with(|| self_access(cx, block)) else {
                return;
            };
            let later = block.stmts[i + 1..]
                .iter()
                .filter_map(stmt_expr)
                .chain(block.expr);
            let mut reentered: Option<Span> = None;
            for e in later {
                match reentered {
                    None => {
                        if refreshes(cx, e, &t, false) {
                            break;
                        }
                        reentered = self.reentry(cx, e, &t, access);
                    }
                    Some(call) => {
                        if refreshes(cx, e, &t, true) {
                            break;
                        }
                        if let Some(at) = reuse_in(cx, e, &t) {
                            let place = dotted(String::from("self"), &t.place);
                            let name = t.name;
                            let msg = match t.fact {
                                Fact::Count | Fact::Flag => format!(
                                    "`{name}` was read from `{place}` before a call that can run other code with access to `self`. It is used here as if `{place}` had not changed"
                                ),
                                Fact::Pointer => format!(
                                    "`{name}` points into `{place}` as it was before a call that can run other code with access to `self`. That code may have moved or freed the buffer by the time `{name}` is used here"
                                ),
                            };
                            emit_with_note(
                                cx,
                                STALE_ACROSS_REENTRY,
                                at,
                                msg,
                                call,
                                format!(
                                    "control leaves this function here, and what runs can reach `{place}`"
                                ),
                                t.fact.help(name, &place),
                            );
                            break;
                        }
                    }
                }
            }
        }
    }
}
