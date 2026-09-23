#![allow(
    clippy::indexing_slicing,
    reason = "graph algorithm indices are always in-bounds"
)]

//! Rewrite support for the `topological_ordering` lint: the target item order
//! of a module and the source text that moves with each item.

use core::cmp::Reverse;
use std::collections::BinaryHeap;

/// Ordering-relevant facts about one module item, indexed like the lint's item
/// list (source order).
pub struct OrderNode {
    /// `const` / `static` items are exempt as ordering targets.
    pub is_const_or_static: bool,
    /// Index of the self type when the item is an impl of a type in the same
    /// module. Such an impl is placed directly after its type.
    pub attached_to: Option<usize>,
}

/// Return the item indices in the order that satisfies the lint.
///
/// `refs` are `(from, to)` edges after impl remapping: `from` must precede
/// `to` unless both share a strongly connected component. Components are
/// emitted in topological order (Kahn's algorithm), preferring the component
/// whose first item comes earliest in the source, so the result stays close
/// to the original order. Items within a component keep their relative order,
/// and each type is followed by its attached impls in their original order.
pub fn desired_order(nodes: &[OrderNode], refs: &[(usize, usize)], item_to_scc: &[usize]) -> Vec<usize> {
    let n = nodes.len();
    let scc_count = item_to_scc.iter().max().map_or(0, |max| max + 1);

    let mut members = vec![Vec::new(); scc_count];
    let mut attached = vec![Vec::new(); n];
    for (idx, node) in nodes.iter().enumerate() {
        match node.attached_to {
            Some(owner) => attached[owner].push(idx),
            None => members[item_to_scc[idx]].push(idx),
        }
    }

    let mut edges: Vec<(usize, usize)> = refs
        .iter()
        .filter(|&&(_, to)| !nodes[to].is_const_or_static)
        .map(|&(from, to)| (item_to_scc[from], item_to_scc[to]))
        .filter(|(from, to)| from != to)
        .collect();
    edges.sort_unstable();
    edges.dedup();

    let mut successors = vec![Vec::new(); scc_count];
    let mut in_degree = vec![0usize; scc_count];
    for &(from, to) in &edges {
        successors[from].push(to);
        in_degree[to] += 1;
    }

    let mut ready: BinaryHeap<Reverse<(usize, usize)>> = (0..scc_count)
        .filter(|&scc| in_degree[scc] == 0)
        .filter_map(|scc| members[scc].first().map(|&first| Reverse((first, scc))))
        .collect();

    let mut order = Vec::with_capacity(n);
    while let Some(Reverse((_, scc))) = ready.pop() {
        for &idx in &members[scc] {
            order.push(idx);
            order.extend_from_slice(&attached[idx]);
        }
        for &next in &successors[scc] {
            in_degree[next] -= 1;
            if in_degree[next] == 0
                && let Some(&first) = members[next].first()
            {
                ready.push(Reverse((first, next)));
            }
        }
    }
    order
}

/// Return the byte offset in `gap` where the text attached to the following
/// item begins, or `None` when `gap` holds anything other than whitespace,
/// comments, and attributes.
///
/// `gap` is the source between the end of the previous item (or the start of
/// the module body) and the start of the item. The attached text is the run of
/// outer attributes and comments directly above the item, not separated from
/// it by a blank line and not preceded by an inner attribute or inner doc
/// comment (those belong to the module). Without such a run, the offset is
/// `gap.len()`.
pub fn attached_start(gap: &str) -> Option<usize> {
    let pieces = scan_gap(gap)?;
    let mut start = gap.len();
    for piece in pieces.iter().rev() {
        if piece.inner || gap[piece.end..start].matches('\n').count() > 1 {
            break;
        }
        start = piece.start;
    }
    Some(start)
}

struct GapPiece {
    start: usize,
    end: usize,
    inner: bool,
}

fn scan_gap(gap: &str) -> Option<Vec<GapPiece>> {
    let bytes = gap.as_bytes();
    let mut pieces = Vec::new();
    let mut pos = 0;
    while pos < bytes.len() {
        if bytes[pos].is_ascii_whitespace() {
            pos += 1;
            continue;
        }
        let rest = &gap[pos..];
        let (end, inner) = if rest.starts_with("//") {
            let end = rest.find('\n').map_or(gap.len(), |newline| pos + newline);
            (end, rest.starts_with("//!"))
        } else if rest.starts_with("/*") {
            (block_comment_end(bytes, pos)?, rest.starts_with("/*!"))
        } else if rest.starts_with("#![") {
            (attribute_end(bytes, pos + 2)?, true)
        } else if rest.starts_with("#[") {
            (attribute_end(bytes, pos + 1)?, false)
        } else {
            return None;
        };
        pieces.push(GapPiece { start: pos, end, inner });
        pos = end;
    }
    Some(pieces)
}

/// Return the offset just past the `*/` that closes the (nestable) block
/// comment opening at `open`.
fn block_comment_end(bytes: &[u8], open: usize) -> Option<usize> {
    let mut depth = 0usize;
    let mut pos = open;
    while pos + 1 < bytes.len() {
        match (bytes[pos], bytes[pos + 1]) {
            (b'/', b'*') => {
                depth += 1;
                pos += 2;
            }
            (b'*', b'/') => {
                depth -= 1;
                pos += 2;
                if depth == 0 {
                    return Some(pos);
                }
            }
            _ => pos += 1,
        }
    }
    None
}

/// Return the offset just past the `]` that closes the attribute bracket at
/// `open`, skipping brackets inside string and character literals.
fn attribute_end(bytes: &[u8], open: usize) -> Option<usize> {
    let mut depth = 0usize;
    let mut pos = open;
    while pos < bytes.len() {
        match bytes[pos] {
            b'[' => depth += 1,
            b']' => {
                depth -= 1;
                if depth == 0 {
                    return Some(pos + 1);
                }
            }
            b'"' => pos = string_end(bytes, pos)?,
            b'r' if !is_ident_byte(bytes[pos - 1]) && matches!(bytes.get(pos + 1), Some(b'"' | b'#')) => {
                pos = raw_string_end(bytes, pos)?;
            }
            b'\'' => pos = char_end(bytes, pos),
            _ => {}
        }
        pos += 1;
    }
    None
}

/// Return the offset of the closing `"` of the string starting at `open`.
fn string_end(bytes: &[u8], open: usize) -> Option<usize> {
    let mut pos = open + 1;
    while pos < bytes.len() {
        match bytes[pos] {
            b'\\' => pos += 2,
            b'"' => return Some(pos),
            _ => pos += 1,
        }
    }
    None
}

/// Return the offset of the last byte of the raw string starting at the `r`
/// at `open`.
fn raw_string_end(bytes: &[u8], open: usize) -> Option<usize> {
    let hashes = bytes[open + 1..].iter().take_while(|&&byte| byte == b'#').count();
    let quote = open + 1 + hashes;
    if bytes.get(quote) != Some(&b'"') {
        return Some(open);
    }
    let mut pos = quote + 1;
    while pos < bytes.len() {
        if bytes[pos] == b'"' && bytes[pos + 1..].iter().take(hashes).filter(|&&byte| byte == b'#').count() == hashes {
            return Some(pos + hashes);
        }
        pos += 1;
    }
    None
}

/// Return the offset of the closing `'` of a character literal at `open`, or
/// `open` itself for a lifetime.
fn char_end(bytes: &[u8], open: usize) -> usize {
    match bytes.get(open + 1) {
        Some(b'\\') => bytes[open + 2..]
            .iter()
            .position(|&byte| byte == b'\'')
            .map_or(open, |offset| open + 2 + offset),
        Some(_) if bytes.get(open + 2) == Some(&b'\'') => open + 2,
        _ => open,
    }
}

fn is_ident_byte(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || byte == b'_'
}

#[cfg(test)]
mod tests {
    use super::{OrderNode, attached_start, desired_order};

    fn node() -> OrderNode {
        OrderNode {
            is_const_or_static: false,
            attached_to: None,
        }
    }

    #[test]
    fn test_attached_start_whitespace_only() {
        let gap = "\n\n";
        assert_eq!(Some(gap.len()), attached_start(gap));
    }

    #[test]
    fn test_attached_start_doc_and_attributes() {
        let gap = "\n\n/// Doc.\n#[derive(Debug)]\n#[cfg_attr(\n    dylint_lib = \"x\",\n    expect(y, reason = \"a ] b\")\n)]\n";
        assert_eq!(Some(2), attached_start(gap));
    }

    #[test]
    fn test_attached_start_stops_at_blank_line() {
        let gap = "\n// SPDX-License-Identifier: MIT\n\n#[inline]\n";
        assert_eq!(Some(34), attached_start(gap));
    }

    #[test]
    fn test_attached_start_stops_at_inner_attribute() {
        let gap = "//! Module doc.\n#![allow(x)]\n#[inline]\n";
        assert_eq!(Some(29), attached_start(gap));
    }

    #[test]
    fn test_attached_start_raw_string_and_block_comment() {
        let gap = "\n/* a /* nested */ b */\n#[doc = r#\"x ] \"# ]\n";
        assert_eq!(Some(1), attached_start(gap));
    }

    #[test]
    fn test_attached_start_rejects_code() {
        assert_eq!(None, attached_start("\n};\n#[inline]\n"));
    }

    #[test]
    fn test_attached_start_rejects_unclosed_attribute() {
        assert_eq!(None, attached_start("\n#[inline\n"));
    }

    #[test]
    fn test_desired_order_moves_callee_below_caller() {
        let nodes = [node(), node(), node()];
        let refs = [(2, 0)];
        let item_to_scc = [0, 1, 2];
        assert_eq!(vec![1, 2, 0], desired_order(&nodes, &refs, &item_to_scc));
    }

    #[test]
    fn test_desired_order_keeps_valid_order() {
        let nodes = [node(), node(), node()];
        let refs = [(0, 1), (1, 2)];
        let item_to_scc = [0, 1, 2];
        assert_eq!(vec![0, 1, 2], desired_order(&nodes, &refs, &item_to_scc));
    }

    #[test]
    fn test_desired_order_ignores_edges_into_const() {
        let nodes = [
            OrderNode {
                is_const_or_static: true,
                attached_to: None,
            },
            node(),
        ];
        let refs = [(1, 0)];
        let item_to_scc = [0, 1];
        assert_eq!(vec![0, 1], desired_order(&nodes, &refs, &item_to_scc));
    }

    #[test]
    fn test_desired_order_places_impl_after_type() {
        let nodes = [
            node(),
            node(),
            OrderNode {
                is_const_or_static: false,
                attached_to: Some(0),
            },
        ];
        let refs = [];
        let item_to_scc = [0, 1, 2];
        assert_eq!(vec![0, 2, 1], desired_order(&nodes, &refs, &item_to_scc));
    }

    #[test]
    fn test_desired_order_keeps_cycle_members_together() {
        let nodes = [node(), node(), node(), node()];
        let refs = [(1, 3), (3, 1), (2, 1)];
        let item_to_scc = [0, 1, 2, 1];
        assert_eq!(vec![0, 2, 1, 3], desired_order(&nodes, &refs, &item_to_scc));
    }
}
