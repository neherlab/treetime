use std::env;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clippy_utils::diagnostics::{span_lint_and_help, span_lint_and_sugg};
use rustc_ast::{AttrKind, Attribute, Crate, Item, ItemKind, MetaItemInner, VariantData};
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
            _ => {},
        }
    }

    fn matched_sources(&self, item: &Item) -> Vec<&RenderSource> {
        let derives = derive_names(&item.attrs);
        self.config
            .rendered
            .iter()
            .filter(|source| {
                source.derives.iter().any(|name| derives.contains(name))
                    || item.attrs.iter().any(|attr| source.attrs.contains(&attr_path(attr)))
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
        let mut offset = strip_shebang(src).unwrap_or(0);
        for token in tokenize(&src[offset..], FrontmatterAllowed::No) {
            let start = offset;
            offset += token.len as usize;
            let doc = match token.kind {
                TokenKind::LineComment { doc_style } | TokenKind::BlockComment { doc_style, .. } => {
                    doc_style.is_some()
                },
                _ => continue,
            };
            let lo = file.start_pos + BytePos(start as u32);
            let hi = file.start_pos + BytePos(offset as u32);
            if doc && self.help_docs.lock().contains(&(lo, hi)) {
                continue;
            }
            if !doc && self.is_allowed_comment(&src[start..offset]) {
                continue;
            }
            report(cx, deletion_span(file, src, start, offset), doc);
        }
    }

    fn is_allowed_comment(&self, text: &str) -> bool {
        let inner = comment_inner_text(text);
        self.config
            .allowed_comment_prefixes
            .iter()
            .any(|prefix| inner.starts_with(prefix.as_str()))
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

fn derive_names(attrs: &[Attribute]) -> Vec<String> {
    attrs
        .iter()
        .filter(|attr| attr.has_name(sym::derive))
        .filter_map(Attribute::meta_item_list)
        .flatten()
        .filter_map(|entry| match entry {
            MetaItemInner::MetaItem(meta) => meta.path.segments.last().map(|segment| segment.ident.to_string()),
            MetaItemInner::Lit(_) => None,
        })
        .collect()
}

fn attr_path(attr: &Attribute) -> String {
    match &attr.kind {
        AttrKind::Normal(normal) => normal
            .item
            .path
            .segments
            .iter()
            .map(|segment| segment.ident.to_string())
            .collect::<Vec<_>>()
            .join("::"),
        AttrKind::DocComment(..) => String::new(),
    }
}

fn sources_rendering<'a>(sources: &[&'a RenderSource], granularity: Granularity) -> Vec<&'a RenderSource> {
    sources
        .iter()
        .copied()
        .filter(|source| source.renders.contains(&granularity))
        .collect()
}

fn is_skipped(source: &RenderSource, attrs: &[Attribute]) -> bool {
    attrs.iter().any(|attr| {
        let path = attr_path(attr);
        source.skip_attrs.iter().any(|spec| {
            let (outer, word) = split_spec(spec);
            path == outer && word.is_none_or(|word| attr_list_has_word(attr, word))
        })
    })
}

fn split_spec(spec: &str) -> (&str, Option<&str>) {
    match spec.split_once('(') {
        Some((outer, rest)) => (outer, Some(rest.trim_end_matches(')'))),
        None => (spec, None),
    }
}

fn attr_list_has_word(attr: &Attribute, word: &str) -> bool {
    attr.meta_item_list().is_some_and(|list| {
        list.iter().any(|entry| match entry {
            MetaItemInner::MetaItem(meta) => meta.is_word() && meta.has_name(Symbol::intern(word)),
            MetaItemInner::Lit(_) => false,
        })
    })
}
