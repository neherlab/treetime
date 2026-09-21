use std::cell::RefCell;
use std::env;
use std::path::{Path, PathBuf};

use clippy_utils::diagnostics::{span_lint_and_help, span_lint_and_sugg};
use rustc_ast::{AttrKind, Attribute, Crate, Item, ItemKind, MetaItemInner, VariantData};
use rustc_data_structures::fx::FxHashSet;
use rustc_errors::Applicability;
use rustc_lexer::{FrontmatterAllowed, TokenKind, strip_shebang, tokenize};
use rustc_lint::{EarlyContext, EarlyLintPass, LintContext as _};
use rustc_span::def_id::LOCAL_CRATE;
use rustc_span::{BytePos, FileName, SourceFile, Span, SyntaxContext, sym};

use crate::config::NoCommentsConfig;

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

thread_local! {
    static HELP_DOCS: RefCell<FxHashSet<(BytePos, BytePos)>> = RefCell::new(FxHashSet::default());
}

pub struct HelpDocCollector {
    config: NoCommentsConfig,
}

impl HelpDocCollector {
    pub fn new() -> Self {
        Self {
            config: dylint_linting::config_or_default("no_comments"),
        }
    }

    fn record_help_docs(&self, cx: &EarlyContext<'_>, item: &Item) {
        let derives = derive_names(&item.attrs);
        let item_help = derives.iter().any(|name| self.config.help_derives.contains(name))
            || item
                .attrs
                .iter()
                .any(|attr| self.config.help_attrs.contains(&attr_path(attr)));
        let variant_help = derives
            .iter()
            .any(|name| self.config.help_variant_derives.contains(name));
        if !item_help && !variant_help {
            return;
        }
        if item_help {
            self.record_doc_block(cx, &item.attrs);
        }
        match &item.kind {
            ItemKind::Struct(_, _, data) => self.record_variant_data(cx, data),
            ItemKind::Enum(_, _, def) => {
                for variant in &def.variants {
                    self.record_doc_block(cx, &variant.attrs);
                    self.record_variant_data(cx, &variant.data);
                }
            },
            _ => {},
        }
    }

    fn record_variant_data(&self, cx: &EarlyContext<'_>, data: &VariantData) {
        let fields = match data {
            VariantData::Struct { fields, .. } | VariantData::Tuple(fields, _) => fields.as_slice(),
            VariantData::Unit(_) => &[],
        };
        for field in fields {
            self.record_doc_block(cx, &field.attrs);
        }
    }

    fn record_doc_block(&self, _cx: &EarlyContext<'_>, attrs: &[Attribute]) {
        let docs = attrs.iter().filter(|attr| attr.doc_str().is_some()).collect::<Vec<_>>();
        let (Some(first), Some(last)) = (docs.first(), docs.last()) else {
            return;
        };
        HELP_DOCS.with(|help| {
            let mut help = help.borrow_mut();
            for attr in &docs {
                help.insert((attr.span.lo(), attr.span.hi()));
            }
        });
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
            _cx,
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
    fn check_item(&mut self, cx: &EarlyContext<'_>, item: &Item) {
        self.record_help_docs(cx, item);
    }
}

pub struct NoComments {
    explicit_docs: Vec<Span>,
    out_dir: Option<PathBuf>,
}

impl NoComments {
    pub fn new() -> Self {
        Self {
            explicit_docs: Vec::new(),
            out_dir: env::var_os("OUT_DIR").map(PathBuf::from),
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
            if doc && HELP_DOCS.with(|help| help.borrow().contains(&(lo, hi))) {
                continue;
            }
            report(cx, deletion_span(file, src, start, offset), doc);
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
                    .is_some_and(|path| !is_generated_path(path, self.out_dir.as_deref())),
                _ => false,
            })
            .cloned()
            .collect::<Vec<_>>();
        for file in &files {
            self.check_file(cx, file);
        }
        for span in &self.explicit_docs {
            let file = source_map.lookup_source_file(span.lo());
            if HELP_DOCS.with(|help| help.borrow().contains(&(span.lo(), span.hi()))) {
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

fn is_generated_path(path: &Path, out_dir: Option<&Path>) -> bool {
    out_dir.is_some_and(|out_dir| path.starts_with(out_dir))
}

fn report(cx: &EarlyContext<'_>, delete: Span, doc: bool) {
    let (msg, help) = if doc {
        (
            "doc comment",
            "the project keeps doc comments only where clap, schemars, or utoipa renders them",
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
