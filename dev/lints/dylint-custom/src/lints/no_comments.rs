use std::env;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clippy_utils::diagnostics::{span_lint_and_help, span_lint_and_sugg};
use rustc_ast::{self as ast, AttrKind, Attribute, Crate, Item, ItemKind, MetaItemInner, VariantData, VisibilityKind};
use rustc_data_structures::fx::FxHashSet;
use rustc_data_structures::sync::Lock;
use rustc_errors::Applicability;
use rustc_lexer::{FrontmatterAllowed, TokenKind, strip_shebang, tokenize};
use rustc_lint::{EarlyContext, EarlyLintPass, LintContext as _};
use rustc_span::def_id::LOCAL_CRATE;
use rustc_span::{BytePos, FileName, SourceFile, Span, Symbol, SyntaxContext, sym};

use crate::config::{Granularity, NoCommentsConfig, RenderSource};

rustc_session::declare_lint! {
    pub NO_COMMENTS,
    Deny,
    "comment in Rust source -- the project keeps explanations in names, types, tests, and the knowledge base"
}

rustc_session::declare_lint! {
    pub DOC_COMMENT_LIMIT,
    Deny,
    "rendered doc comment over the size limit -- keep the summary, move the rest to the knowledge base"
}

pub type HelpDocs = Arc<Lock<FxHashSet<(BytePos, BytePos)>>>;

pub fn new_help_docs() -> HelpDocs {
    Arc::new(Lock::new(FxHashSet::default()))
}

pub struct HelpDocCollector {
    config: NoCommentsConfig,
    help_docs: HelpDocs,
}

impl HelpDocCollector {
    pub fn new(help_docs: HelpDocs) -> Self {
        Self {
            config: dylint_linting::config_or_default("no_comments"),
            help_docs,
        }
    }

    fn record_help_docs(&self, cx: &EarlyContext<'_>, item: &Item) {
        let sources = self.matched_sources(item);
        if sources.is_empty() {
            return;
        }
        if sources.iter().any(|source| source.renders.contains(&Granularity::Item)) {
            self.record_doc_block(cx, &item.attrs);
        }
        match &item.kind {
            ItemKind::Struct(_, _, data) => {
                let rendering = sources_rendering(&sources, Granularity::Fields);
                self.record_field_docs(cx, &rendering, data);
            },
            ItemKind::Enum(_, _, def) => {
                for variant in &def.variants {
                    let rendering = sources_rendering(&sources, Granularity::Variants)
                        .into_iter()
                        .filter(|source| !is_skipped(source, &variant.attrs))
                        .collect::<Vec<_>>();
                    if rendering.is_empty() {
                        continue;
                    }
                    self.record_doc_block(cx, &variant.attrs);
                    self.record_field_docs(cx, &rendering, &variant.data);
                }
            },
            ItemKind::Impl(imp) => {
                let rendering = sources_rendering(&sources, Granularity::Methods);
                for assoc in &imp.items {
                    let shown = matches!(assoc.vis.kind, VisibilityKind::Public)
                        && rendering.iter().any(|source| !is_skipped(source, &assoc.attrs));
                    if shown {
                        self.record_doc_block(cx, &assoc.attrs);
                    }
                }
            },
            _ => {},
        }
    }

    fn matched_sources(&self, item: &Item) -> Vec<&RenderSource> {
        let attrs = item.attrs.iter().flat_map(attr_views).collect::<Vec<_>>();
        let derives = derive_names(&attrs);
        self.config
            .rendered
            .iter()
            .filter(|source| {
                source.derives.iter().any(|name| derives.contains(name))
                    || attrs.iter().any(|view| source.attrs.contains(&view.path))
            })
            .collect()
    }

    fn record_field_docs(&self, cx: &EarlyContext<'_>, rendering: &[&RenderSource], data: &VariantData) {
        if rendering.is_empty() {
            return;
        }
        let fields = match data {
            VariantData::Struct { fields, .. } | VariantData::Tuple(fields, _) => fields.as_slice(),
            VariantData::Unit(_) => &[],
        };
        for field in fields {
            if rendering.iter().any(|source| !is_skipped(source, &field.attrs)) {
                self.record_doc_block(cx, &field.attrs);
            }
        }
    }

    fn record_doc_block(&self, cx: &EarlyContext<'_>, attrs: &[Attribute]) {
        let docs = attrs.iter().filter(|attr| attr.doc_str().is_some()).collect::<Vec<_>>();
        let (Some(first), Some(last)) = (docs.first(), docs.last()) else {
            return;
        };
        {
            let mut help = self.help_docs.lock();
            for attr in &docs {
                help.insert((attr.span.lo(), attr.span.hi()));
            }
        }
        let lines = docs
            .iter()
            .filter_map(|attr| attr.doc_str())
            .flat_map(|text| {
                text.to_string()
                    .lines()
                    .map(str::trim)
                    .map(str::to_owned)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let text = lines.join("\n");
        let paragraphs = text
            .split("\n\n")
            .filter(|paragraph| !paragraph.trim().is_empty())
            .count();
        let over = [
            (
                text.chars().count() > self.config.max_chars,
                format!("{} characters", self.config.max_chars),
            ),
            (
                lines.len() > self.config.max_lines,
                format!("{} lines", self.config.max_lines),
            ),
            (
                paragraphs > self.config.max_paragraphs,
                format!("{} paragraphs", self.config.max_paragraphs),
            ),
        ]
        .into_iter()
        .filter_map(|(over, limit)| over.then_some(limit))
        .collect::<Vec<_>>();
        if over.is_empty() {
            return;
        }
        span_lint_and_help(
            cx,
            DOC_COMMENT_LIMIT,
            first.span.to(last.span),
            format!("rendered doc comment is over {}", over.join(", ")),
            None,
            "keep the summary here and move the rest to the knowledge base",
        );
    }
}

rustc_session::impl_lint_pass!(HelpDocCollector => [DOC_COMMENT_LIMIT]);

impl EarlyLintPass for HelpDocCollector {
    fn check_crate(&mut self, _cx: &EarlyContext<'_>, _krate: &Crate) {
        self.help_docs.lock().clear();
    }

    fn check_item(&mut self, cx: &EarlyContext<'_>, item: &Item) {
        self.record_help_docs(cx, item);
    }
}

pub struct NoComments {
    config: NoCommentsConfig,
    explicit_docs: Vec<Span>,
    cargo_out_dir: Option<PathBuf>,
    help_docs: HelpDocs,
}

impl NoComments {
    pub fn new(help_docs: HelpDocs) -> Self {
        Self {
            config: dylint_linting::config_or_default("no_comments"),
            explicit_docs: Vec::new(),
            cargo_out_dir: env::var_os("OUT_DIR").map(PathBuf::from),
            help_docs,
        }
    }

    fn check_file(&self, cx: &EarlyContext<'_>, file: &SourceFile) {
        let Some(src) = file.src.as_deref() else {
            return;
        };
        for comment in scan_comments(src, &self.config.allowed_comment_prefixes) {
            if comment.kept {
                continue;
            }
            let lo = file.start_pos + BytePos(comment.start as u32);
            let hi = file.start_pos + BytePos(comment.end as u32);
            if comment.doc && self.help_docs.lock().contains(&(lo, hi)) {
                continue;
            }
            report(cx, deletion_span(file, src, comment.start, comment.end), comment.doc);
        }
    }
}

rustc_session::impl_lint_pass!(NoComments => [NO_COMMENTS]);

impl EarlyLintPass for NoComments {
    fn check_crate(&mut self, _cx: &EarlyContext<'_>, _krate: &Crate) {
        self.explicit_docs.clear();
    }

    fn check_attribute(&mut self, _cx: &EarlyContext<'_>, attr: &Attribute) {
        if attr.has_name(sym::doc)
            && !attr.is_doc_comment()
            && attr.value_str().is_some()
            && !attr.span.from_expansion()
        {
            self.explicit_docs.push(attr.span);
        }
    }

    fn check_crate_post(&mut self, cx: &EarlyContext<'_>, _krate: &Crate) {
        let source_map = cx.sess().source_map();
        let files = source_map
            .files()
            .iter()
            .filter(|file| file.cnum == LOCAL_CRATE)
            .filter(|file| match &file.name {
                FileName::Real(name) => name
                    .local_path()
                    .is_some_and(|path| !is_generated_path(path, self.cargo_out_dir.as_deref())),
                _ => false,
            })
            .cloned()
            .collect::<Vec<_>>();
        for file in &files {
            self.check_file(cx, file);
        }
        for span in &self.explicit_docs {
            let file = source_map.lookup_source_file(span.lo());
            if self.help_docs.lock().contains(&(span.lo(), span.hi())) {
                continue;
            }
            let Some(src) = file.src.as_deref() else {
                continue;
            };
            let start = (span.lo() - file.start_pos).0 as usize;
            let end = (span.hi() - file.start_pos).0 as usize;
            report(cx, deletion_span(&file, src, start, end), true);
        }
    }
}

fn scan_comments(src: &str, allowed_prefixes: &[String]) -> Vec<ScannedComment> {
    let mut offset = strip_shebang(src).unwrap_or(0);
    let mut kept_line_end = None;
    let mut comments = Vec::new();
    for token in tokenize(&src[offset..], FrontmatterAllowed::No) {
        let start = offset;
        offset += token.len as usize;
        let (doc, line) = match token.kind {
            TokenKind::LineComment { doc_style } => (doc_style.is_some(), true),
            TokenKind::BlockComment { doc_style, .. } => (doc_style.is_some(), false),
            TokenKind::Whitespace => continue,
            _ => {
                kept_line_end = None;
                continue;
            },
        };
        let continues_kept =
            kept_line_end.is_some_and(|end: usize| src[end..start].matches('\n').count() == 1);
        let kept = !doc && (continues_kept || has_allowed_prefix(&src[start..offset], allowed_prefixes));
        kept_line_end = (kept && line).then_some(offset);
        comments.push(ScannedComment {
            start,
            end: offset,
            doc,
            kept,
        });
    }
    comments
}

fn has_allowed_prefix(text: &str, allowed_prefixes: &[String]) -> bool {
    let inner = comment_inner_text(text);
    allowed_prefixes.iter().any(|prefix| inner.starts_with(prefix.as_str()))
}

struct ScannedComment {
    start: usize,
    end: usize,
    doc: bool,
    kept: bool,
}

fn is_generated_path(path: &Path, cargo_out_dir: Option<&Path>) -> bool {
    cargo_out_dir.is_some_and(|cargo_out_dir| path.starts_with(cargo_out_dir))
}

fn report(cx: &EarlyContext<'_>, delete: Span, doc: bool) {
    let (msg, help) = if doc {
        (
            "doc comment",
            "keep a doc comment only where a configured render source uses it",
        )
    } else {
        (
            "comment",
            "the project keeps no comments; put the fact in a name, a type, a test, or the knowledge base",
        )
    };
    span_lint_and_sugg(
        cx,
        NO_COMMENTS,
        delete,
        msg,
        help,
        String::new(),
        Applicability::MachineApplicable,
    );
}

fn deletion_span(file: &SourceFile, src: &str, start: usize, end: usize) -> Span {
    let (delete_start, delete_end) = deletion_range(src, start, end);
    Span::new(
        file.start_pos + BytePos(delete_start as u32),
        file.start_pos + BytePos(delete_end as u32),
        SyntaxContext::root(),
        None,
    )
}

fn deletion_range(src: &str, start: usize, end: usize) -> (usize, usize) {
    let line_start = src[..start].rfind('\n').map_or(0, |pos| pos + 1);
    if src[line_start..start].trim().is_empty() {
        let rest = &src[end..];
        let newline = if rest.starts_with("\r\n") {
            2
        } else {
            usize::from(rest.starts_with('\n'))
        };
        (line_start, end + newline)
    } else {
        (src[..start].trim_end().len(), end)
    }
}

fn comment_inner_text(text: &str) -> &str {
    let text = text.trim();
    let text = text
        .strip_prefix("/*")
        .map_or(text, |rest| rest.trim_end_matches("*/"));
    text.trim_start_matches('/').trim()
}

/// An attribute as it applies after `cfg_attr` expansion: its path and its
/// argument list. The pre-expansion pass sees `#[cfg_attr(pred, attr, ...)]`
/// unexpanded, so the attributes inside it are taken as present whatever the
/// predicate, because a doc comment that any configuration renders is
/// user-facing.
struct AttrView {
    path: String,
    list: Vec<MetaItemInner>,
}

fn attr_views(attr: &Attribute) -> Vec<AttrView> {
    if attr.has_name(sym::cfg_attr) {
        return attr
            .meta_item_list()
            .into_iter()
            .flatten()
            .skip(1)
            .filter_map(|entry| match entry {
                MetaItemInner::MetaItem(meta) => Some(AttrView {
                    path: path_string(&meta.path),
                    list: meta.meta_item_list().map(<[_]>::to_vec).unwrap_or_default(),
                }),
                MetaItemInner::Lit(_) => None,
            })
            .collect();
    }
    match &attr.kind {
        AttrKind::Normal(normal) => vec![AttrView {
            path: path_string(&normal.item.path),
            list: attr.meta_item_list().map(Vec::from).unwrap_or_default(),
        }],
        AttrKind::DocComment(..) => vec![],
    }
}

fn path_string(path: &ast::Path) -> String {
    path.segments
        .iter()
        .map(|segment| segment.ident.to_string())
        .collect::<Vec<_>>()
        .join("::")
}

fn derive_names(attrs: &[AttrView]) -> Vec<String> {
    attrs
        .iter()
        .filter(|view| view.path == "derive")
        .flat_map(|view| &view.list)
        .filter_map(|entry| match entry {
            MetaItemInner::MetaItem(meta) => meta.path.segments.last().map(|segment| segment.ident.to_string()),
            MetaItemInner::Lit(_) => None,
        })
        .collect()
}

fn sources_rendering<'a>(sources: &[&'a RenderSource], granularity: Granularity) -> Vec<&'a RenderSource> {
    sources
        .iter()
        .copied()
        .filter(|source| source.renders.contains(&granularity))
        .collect()
}

fn is_skipped(source: &RenderSource, attrs: &[Attribute]) -> bool {
    attrs.iter().flat_map(attr_views).any(|view| {
        source.skip_attrs.iter().any(|spec| {
            let (outer, word) = split_spec(spec);
            view.path == outer && word.is_none_or(|word| list_has_word(&view.list, word))
        })
    })
}

fn split_spec(spec: &str) -> (&str, Option<&str>) {
    match spec.split_once('(') {
        Some((outer, rest)) => (outer, Some(rest.trim_end_matches(')'))),
        None => (spec, None),
    }
}

fn list_has_word(list: &[MetaItemInner], word: &str) -> bool {
    list.iter().any(|entry| match entry {
        MetaItemInner::MetaItem(meta) => meta.is_word() && meta.has_name(Symbol::intern(word)),
        MetaItemInner::Lit(_) => false,
    })
}

#[cfg(test)]
mod tests {
    use super::scan_comments;

    fn kept(src: &str) -> Vec<(String, bool)> {
        let prefixes = ["SAFETY:", "TODO", "HACK"].map(str::to_owned);
        scan_comments(src, &prefixes)
            .into_iter()
            .map(|comment| (src[comment.start..comment.end].to_owned(), comment.kept))
            .collect()
    }

    #[test]
    fn test_scan_comments_keeps_marker_line() {
        assert_eq!(
            vec![("// SAFETY: ascii only".to_owned(), true)],
            kept("let x = 1;\n// SAFETY: ascii only\nunsafe {}\n")
        );
    }

    #[test]
    fn test_scan_comments_rejects_plain_line() {
        assert_eq!(vec![("// explains x".to_owned(), false)], kept("// explains x\nlet x = 1;\n"));
    }

    #[test]
    fn test_scan_comments_keeps_continuation_lines() {
        let src = "  // TODO: enable when fixed\n  // #[case::a(\"a\")] // slow\n  // #[case::b(\"b\")]\nfn f() {}\n";
        assert_eq!(
            vec![
                ("// TODO: enable when fixed".to_owned(), true),
                ("// #[case::a(\"a\")] // slow".to_owned(), true),
                ("// #[case::b(\"b\")]".to_owned(), true),
            ],
            kept(src)
        );
    }

    #[test]
    fn test_scan_comments_blank_line_ends_group() {
        assert_eq!(
            vec![("// HACK: guess".to_owned(), true), ("// unrelated".to_owned(), false)],
            kept("// HACK: guess\n\n// unrelated\n")
        );
    }

    #[test]
    fn test_scan_comments_code_ends_group() {
        assert_eq!(
            vec![("// TODO: x".to_owned(), true), ("// unrelated".to_owned(), false)],
            kept("// TODO: x\nlet y = 2;\n// unrelated\n")
        );
    }

    #[test]
    fn test_scan_comments_trailing_marker_after_code() {
        assert_eq!(vec![("// TODO: slow".to_owned(), true)], kept("let z = 3; // TODO: slow\n"));
    }

    #[test]
    fn test_scan_comments_marker_block_comment() {
        assert_eq!(vec![("/* TODO: x\n y */".to_owned(), true)], kept("/* TODO: x\n y */\n"));
    }

    #[test]
    fn test_scan_comments_doc_comment_not_kept() {
        assert_eq!(
            vec![("/// TODO: x".to_owned(), false), ("// next".to_owned(), false)],
            kept("/// TODO: x\n// next\nfn f() {}\n")
        );
    }

    #[test]
    fn test_scan_comments_doc_after_marker_not_kept() {
        assert_eq!(
            vec![("// TODO: x".to_owned(), true), ("/// doc".to_owned(), false)],
            kept("// TODO: x\n/// doc\nfn f() {}\n")
        );
    }
}
