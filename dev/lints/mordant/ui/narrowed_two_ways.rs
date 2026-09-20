// A place range-checked at one narrowing must not be truncated with `as` at another.

use std::convert::{TryFrom, TryInto};

struct Buf {
    len: usize,
    offset: u64,
    count: u32,
    hash: u64,
    code: i32,
    small: u16,
    point: i32,
    start: usize,
}

// Flagged: `header_len` checks `len` into u32, this truncates it.
fn wire_len(b: &Buf) -> u32 {
    b.len as u32
}

fn header_len(b: &Buf) -> u32 {
    u32::try_from(b.len).expect("length fits the header")
}

// Flagged: same-width sign flip; `seek_to` checks it, this reinterprets it.
fn tell(b: &Buf) -> i64 {
    b.offset as i64
}

fn seek_to(b: &Buf) -> i64 {
    i64::try_from(b.offset).expect("offset fits off_t")
}

// Flagged: `try_into` is a check like `try_from`.
fn count_byte(b: &Buf) -> u8 {
    b.count as u8
}

fn count_checked(b: &Buf) -> u8 {
    let c: u8 = b.count.try_into().unwrap_or(u8::MAX);
    c
}

// Flagged: one local, checked and truncated in the same body.
fn local_both_ways(n: u64) -> (u16, u16) {
    let checked = u16::try_from(n).unwrap_or(u16::MAX);
    (n as u16, checked)
}

// Fine: the comparison against a constant is the check for the cast after it.
fn compared_then_cast(b: &Buf) -> u32 {
    if b.len <= u32::MAX as usize {
        b.len as u32
    } else {
        u32::MAX
    }
}

// Fine: a range pattern is a check too.
fn matched_then_cast(n: u64) -> u8 {
    let wide = u8::try_from(n).is_err();
    match n {
        0..=255 if !wide => n as u8,
        _ => 0,
    }
}

// Fine: an `if let` range pattern is the same check.
fn if_let_then_cast(n: u64) -> u8 {
    let narrow = u8::try_from(n).ok();
    if let 0..=255 = n {
        n as u8
    } else {
        narrow.unwrap_or(0)
    }
}

// Fine: out through `usize` and back is `u32` again; the operand's type is not
// the place's.
fn round_trip(n: u32) -> (u32, bool) {
    (n as usize as u32, i32::try_from(n).is_err())
}

// Fine: `repr(C)` fixes the field's width; it cannot be redeclared narrow.
#[repr(C)]
struct Header {
    size: u64,
}

fn header_size(h: &Header) -> u32 {
    h.size as u32
}

fn header_size_checked(h: &Header) -> u32 {
    u32::try_from(h.size).unwrap_or(u32::MAX)
}

// Fine: widening loses nothing.
fn widen(b: &Buf) -> u64 {
    b.small as u64 + b.len as u64
}

// Fine: nobody checks `hash`; truncating it everywhere is the convention.
fn bucket(b: &Buf) -> u32 {
    b.hash as u32
}

fn bucket_again(b: &Buf) -> u16 {
    b.hash as u16
}

// Fine: a computed operand is not the place.
fn low_byte(b: &Buf) -> u8 {
    (b.len & 0xff) as u8
}

// Fine: the check into `u8` sits where `code` is already known small; it says
// nothing about whether `code` fits the wider `u32` elsewhere.
fn flag_index(b: &Buf) -> u8 {
    if b.code >= 0x61 && b.code <= 0x7a {
        u8::try_from(b.code).expect("a lowercase letter") - 0x61
    } else {
        0
    }
}

fn shown(b: &Buf) -> char {
    char::from_u32(b.code as u32).unwrap_or('?')
}

// Fine: comparing `point` against a constant excuses the cast beside it, but
// is no evidence against the cast in `is_letter`: what it guards is not said.
fn in_bmp(b: &Buf) -> bool {
    b.point <= 0xFFFF && char::from_u32(b.point as u32).is_some()
}

fn is_letter(b: &Buf) -> bool {
    char::from_u32(b.point as u32).map_or(false, char::is_alphabetic)
}

// Fine: the mask right after the cast keeps the low bits on purpose.
fn bit_in_word(b: &Buf) -> u32 {
    (b.start as u32) & (usize::BITS - 1)
}

fn word_checked(b: &Buf) -> u32 {
    u32::try_from(b.start).expect("small set") / usize::BITS
}

// Fine: a local is a place only within its own body.
fn checks_its_own(n: u64) -> u32 {
    u32::try_from(n).unwrap_or(u32::MAX)
}

fn truncates_its_own(n: u64) -> u32 {
    n as u32
}

fn main() {
    let b = Buf {
        len: 1,
        offset: 2,
        count: 3,
        hash: 4,
        code: 5,
        small: 6,
        point: 11,
        start: 12,
    };
    let _ = wire_len(&b);
    let _ = header_len(&b);
    let _ = tell(&b);
    let _ = seek_to(&b);
    let _ = count_byte(&b);
    let _ = count_checked(&b);
    let _ = local_both_ways(7);
    let _ = compared_then_cast(&b);
    let _ = matched_then_cast(8);
    let _ = if_let_then_cast(8);
    let _ = round_trip(3);
    let h = Header { size: 13 };
    let _ = header_size(&h);
    let _ = header_size_checked(&h);
    let _ = widen(&b);
    let _ = bucket(&b);
    let _ = bucket_again(&b);
    let _ = low_byte(&b);
    let _ = flag_index(&b);
    let _ = shown(&b);
    let _ = in_bmp(&b);
    let _ = is_letter(&b);
    let _ = bit_in_word(&b);
    let _ = word_checked(&b);
    let _ = checks_its_own(9);
    let _ = truncates_its_own(10);
}
