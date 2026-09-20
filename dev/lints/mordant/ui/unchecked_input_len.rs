// A parameter, or a field reached from one, that one path bounds (an ordering
// comparison, possibly of arithmetic on it) and another path, or the same path
// earlier, hands to something that turns it into memory is flagged at the use;
// an assertion about something else in front of the use is not the check. Not
// flagged: uses every path checks first (debug_assert!, `.get(i)?`, a call or
// assert_eq! given `&value` included), arithmetic on the value, a copy the body
// advances or re-assigns, `&value` itself, `self`'s fields, values only tested
// for equality or against a range, never compared, computed, written to, slicing.

pub struct Header {
    pub kind: u8,
    pub len: usize,
    pub count: u32,
}

pub struct Frame {
    pub header: Header,
    pub payload: Vec<u8>,
}

fn sibling_arm_is_flagged<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 {
        if hdr.len > buf.len() {
            return &[];
        }
        return unsafe { buf.get_unchecked(..hdr.len) };
    }
    unsafe { buf.get_unchecked(..hdr.len) }
}

fn use_before_check_is_flagged(len: usize, buf: &[u8]) -> Option<&[u8]> {
    let (head, _) = buf.split_at(len);
    if len > buf.len() {
        return None;
    }
    Some(head)
}

fn nested_field_through_local_is_flagged(f: &Frame) -> Vec<u8> {
    let n = f.header.len;
    let out = Vec::with_capacity(n);
    if n > 1 << 20 {
        return Vec::new();
    }
    out
}

fn conditional_check_then_unconditional_use_is_flagged(hdr: &Header, out: &mut Vec<u8>) {
    if hdr.kind == 2 && hdr.count > 1024 {
        return;
    }
    out.reserve(hdr.count as usize);
}

fn capture_is_flagged(len: usize, buf: &[u8]) -> impl Fn(bool) -> usize + '_ {
    move |strict| {
        if strict && len > buf.len() {
            return 0;
        }
        unsafe { buf.get_unchecked(..len) }.len()
    }
}

fn check_through_arithmetic_is_flagged<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 {
        if hdr.len + 4 > buf.len() {
            return &[];
        }
        return &buf[4..];
    }
    unsafe { buf.get_unchecked(..hdr.len) }
}

fn pointer_offset_is_flagged(off: usize, base: *const u8, cap: usize) -> *const u8 {
    let p = unsafe { base.add(off) };
    if off >= cap {
        return base;
    }
    p
}

fn check_in_one_arm_beside_a_panicking_arm_is_flagged(hdr: &Header, buf: &[u8]) -> usize {
    match hdr.kind {
        0 => {
            if hdr.len > buf.len() {
                return 0;
            }
        }
        1 => {}
        _ => unreachable!(),
    }
    buf.split_at(hdr.len).0.len()
}

fn assertion_about_something_else_before_the_use_is_flagged(
    ok: bool,
    len: usize,
    buf: &[u8],
) -> Option<&[u8]> {
    assert!(ok);
    let (head, _) = buf.split_at(len);
    if len > buf.len() {
        return None;
    }
    Some(head)
}

fn let_else_before_the_use_is_flagged(first: Option<u8>, len: usize, buf: &[u8]) -> Option<&[u8]> {
    let Some(_first) = first else { panic!() };
    let (head, _) = buf.split_at(len);
    if len > buf.len() {
        return None;
    }
    Some(head)
}

fn assertion_beside_the_use_in_an_arm_is_flagged(
    strict: bool,
    len: usize,
    buf: &[u8],
) -> *const u8 {
    if strict {
        assert!(!buf.is_empty());
        let p = unsafe { buf.as_ptr().add(len) };
        if len > buf.len() {
            return buf.as_ptr();
        }
        return p;
    }
    buf.as_ptr()
}

fn judging_call_that_does_not_dominate_is_flagged<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 && (hdr.len > buf.len() || !fits(&hdr.len, buf)) {
        return &[];
    }
    buf.split_at(hdr.len).0
}

unsafe fn read_through_a_copied_pointer_is_flagged(p: *const Header, buf: &[u8]) -> usize {
    let q = p;
    if unsafe { (*q).kind } == 1 && unsafe { (*q).len } > buf.len() {
        return 0;
    }
    buf.split_at(unsafe { (*q).len }).0.len()
}

fn loop_over_a_range_of_it_is_not_a_check_so_is_flagged(len: usize, buf: &[u8]) -> usize {
    let mut sum = 0;
    for i in 0..len {
        sum += i;
    }
    if sum == 0 && len > buf.len() {
        return 0;
    }
    buf.split_at(len).0.len()
}

fn fits(len: &usize, buf: &[u8]) -> bool {
    *len <= buf.len()
}

fn every_path_checked_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.len > buf.len() {
        return &[];
    }
    if hdr.kind == 1 {
        return buf.split_at(hdr.len).0;
    }
    buf.split_at(hdr.len).1
}

fn clamp_dominates_is_fine(len: usize, buf: &[u8]) -> &[u8] {
    let mut n = len;
    if n > buf.len() {
        n = buf.len();
    }
    buf.split_at(n).0
}

fn assert_dominates_is_fine(len: usize, buf: &[u8]) -> &[u8] {
    assert!(len <= buf.len());
    if len == 0 {
        return &[];
    }
    buf.split_at(len).0
}

fn debug_assert_dominates_is_fine(len: usize, buf: &[u8]) -> &[u8] {
    debug_assert!(len <= buf.len());
    if len == 0 {
        return &[];
    }
    buf.split_at(len).0
}

fn conjunction_in_assertion_is_fine<'a>(a: usize, b: usize, buf: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    debug_assert!(a <= buf.len() && b <= buf.len());
    (buf.split_at(a).0, buf.split_at(b).0)
}

fn match_arms_are_fine(len: usize, buf: &[u8]) -> &[u8] {
    match len {
        0 => &[],
        1..=8 => buf.split_at(len).0,
        n => buf.split_at(n.min(buf.len())).0,
    }
}

fn never_compared_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    buf.split_at(hdr.len).0
}

fn through_a_call_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 && hdr.len > buf.len() {
        return &[];
    }
    let n = hdr.len.min(buf.len());
    buf.split_at(n).0
}

fn own_length_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    let n = buf.len();
    if hdr.kind == 1 && n < 4 {
        return &[];
    }
    buf.split_at(n).0
}

fn payload_of_call_is_fine(input: Option<usize>, buf: &[u8]) -> Option<&[u8]> {
    let n = input?;
    let tail = unsafe { buf.get_unchecked(..n) };
    if n > buf.len() {
        return None;
    }
    Some(tail)
}

fn written_place_is_fine<'a>(hdr: &mut Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 && hdr.len > buf.len() {
        hdr.len = buf.len();
    }
    buf.split_at(hdr.len).0
}

fn sibling_field_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 && hdr.count > 8 {
        return &[];
    }
    buf.split_at(hdr.len).0
}

impl Header {
    fn own_field_is_fine<'a>(&self, buf: &'a [u8]) -> &'a [u8] {
        let out = buf.split_at(self.len).0;
        if self.len < 16 {
            return &[];
        }
        out
    }
}

fn equality_is_not_a_bound_so_is_fine(len: usize, buf: &[u8]) -> &[u8] {
    let head = buf.split_at(len).0;
    if len == 0 || len != buf.len() {
        return &[];
    }
    head
}

fn equality_check_still_covers_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 {
        if hdr.len != 8 {
            return &[];
        }
        return buf.split_at(hdr.len).0;
    }
    if hdr.len > buf.len() {
        return &[];
    }
    buf.split_at(hdr.len).0
}

fn judging_call_covers_is_fine(p: usize, buf: &[u8]) -> Option<&[u8]> {
    let first = *buf.get(p)?;
    if first == b'-' {
        return Some(buf.split_at(p).0);
    }
    if p > 4 {
        return None;
    }
    Some(buf.split_at(p).1)
}

fn derived_quantity_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    let end = 4 + hdr.len;
    if hdr.kind == 1 && hdr.len > buf.len() {
        return &[];
    }
    buf.split_at(end).0
}

fn copy_that_is_later_advanced_is_fine(p: usize, buf: &[u8]) -> usize {
    let mut q = p;
    let (head, _) = buf.split_at(q);
    while q < buf.len() && buf[q] != b',' {
        q += 1;
    }
    q + head.len()
}

fn safe_indexing_is_not_a_sink_so_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> (&'a [u8], u8) {
    if hdr.kind == 1 && hdr.len > buf.len() {
        return (&[], 0);
    }
    (&buf[..hdr.len], buf[hdr.count as usize])
}

fn copy_reassigned_on_the_other_path_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    let mut n = hdr.len;
    if hdr.kind == 1 {
        if hdr.len > buf.len() {
            return &[];
        }
    } else {
        n = buf.len();
    }
    buf.split_at(n).0
}

fn copy_clamped_on_the_other_path_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    let mut n = hdr.len;
    if hdr.kind == 1 {
        if n > buf.len() {
            return &[];
        }
    } else {
        n = n.min(buf.len());
    }
    unsafe { buf.get_unchecked(..n) }
}

fn call_given_a_reference_covers_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 && hdr.len > buf.len() {
        return &[];
    }
    if !fits(&hdr.len, buf) {
        return &[];
    }
    buf.split_at(hdr.len).0
}

fn assert_eq_covers_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 && hdr.len > buf.len() {
        return &[];
    }
    assert_eq!(hdr.len, buf.len());
    buf.split_at(hdr.len).0
}

fn cmp_covers_is_fine<'a>(hdr: &Header, buf: &'a [u8]) -> &'a [u8] {
    if hdr.kind == 1 && hdr.len > buf.len() {
        return &[];
    }
    if let std::cmp::Ordering::Greater = hdr.len.cmp(&buf.len()) {
        return &[];
    }
    buf.split_at(hdr.len).0
}

fn address_of_the_value_is_not_the_value_so_is_fine<'a>(
    hdr: &'a Header,
    buf: &[u8],
) -> &'a [usize] {
    if hdr.kind == 1 && hdr.len > buf.len() {
        return &[];
    }
    unsafe { std::slice::from_raw_parts(&hdr.len, 1) }
}

fn assertion_before_a_debug_assert_is_fine(ok: bool, len: usize, buf: &[u8]) -> &[u8] {
    assert!(ok);
    debug_assert!(len <= buf.len());
    if len == 0 {
        return &[];
    }
    buf.split_at(len).0
}

fn range_pattern_is_a_dispatch_so_is_fine(len: usize, buf: &[u8]) -> usize {
    let mut weight = 0;
    match len {
        0 => weight += 1,
        1..=8 => weight += 2,
        _ => weight += 3,
    }
    buf.split_at(len).0.len() + weight
}

unsafe fn write_through_a_copied_pointer_is_fine(p: *mut Header, buf: &[u8]) -> usize {
    let q = p;
    unsafe { (*q).len = buf.len() };
    if unsafe { (*p).kind } == 1 && unsafe { (*p).len } > buf.len() {
        return 0;
    }
    buf.split_at(unsafe { (*p).len }).0.len()
}

fn main() {
    let hdr = Header {
        kind: 1,
        len: 2,
        count: 1,
    };
    let frame = Frame {
        header: Header {
            kind: 0,
            len: 0,
            count: 0,
        },
        payload: Vec::new(),
    };
    let buf = [0u8; 8];
    let mut owned = Header {
        kind: 1,
        len: 99,
        count: 0,
    };
    let whole = Header {
        kind: 0,
        len: buf.len(),
        count: 0,
    };
    let mut out = Vec::new();
    let _ = sibling_arm_is_flagged(&hdr, &buf);
    let _ = use_before_check_is_flagged(2, &buf);
    let _ = nested_field_through_local_is_flagged(&frame);
    let _ = frame.payload.len();
    conditional_check_then_unconditional_use_is_flagged(&hdr, &mut out);
    let _ = capture_is_flagged(2, &buf)(true);
    let _ = check_through_arithmetic_is_flagged(&hdr, &buf);
    let _ = pointer_offset_is_flagged(1, buf.as_ptr(), buf.len());
    let _ = check_in_one_arm_beside_a_panicking_arm_is_flagged(&hdr, &buf);
    let _ = assertion_about_something_else_before_the_use_is_flagged(true, 2, &buf);
    let _ = let_else_before_the_use_is_flagged(Some(0), 2, &buf);
    let _ = assertion_beside_the_use_in_an_arm_is_flagged(true, 2, &buf);
    let _ = judging_call_that_does_not_dominate_is_flagged(&hdr, &buf);
    let _ = unsafe { read_through_a_copied_pointer_is_flagged(&hdr, &buf) };
    let _ = loop_over_a_range_of_it_is_not_a_check_so_is_flagged(2, &buf);
    let _ = every_path_checked_is_fine(&hdr, &buf);
    let _ = clamp_dominates_is_fine(2, &buf);
    let _ = assert_dominates_is_fine(2, &buf);
    let _ = debug_assert_dominates_is_fine(2, &buf);
    let _ = conjunction_in_assertion_is_fine(1, 2, &buf);
    let _ = match_arms_are_fine(2, &buf);
    let _ = never_compared_is_fine(&hdr, &buf);
    let _ = through_a_call_is_fine(&hdr, &buf);
    let _ = own_length_is_fine(&hdr, &buf);
    let _ = payload_of_call_is_fine(Some(2), &buf);
    let _ = written_place_is_fine(&mut owned, &buf);
    let _ = sibling_field_is_fine(&hdr, &buf);
    let _ = hdr.own_field_is_fine(&buf);
    let _ = equality_is_not_a_bound_so_is_fine(2, &buf);
    let _ = equality_check_still_covers_is_fine(&hdr, &buf);
    let _ = judging_call_covers_is_fine(1, &buf);
    let _ = derived_quantity_is_fine(&hdr, &buf);
    let _ = copy_that_is_later_advanced_is_fine(1, &buf);
    let _ = safe_indexing_is_not_a_sink_so_is_fine(&hdr, &buf);
    let _ = copy_reassigned_on_the_other_path_is_fine(&hdr, &buf);
    let _ = copy_clamped_on_the_other_path_is_fine(&hdr, &buf);
    let _ = call_given_a_reference_covers_is_fine(&hdr, &buf);
    let _ = assert_eq_covers_is_fine(&whole, &buf);
    let _ = cmp_covers_is_fine(&hdr, &buf);
    let _ = address_of_the_value_is_not_the_value_so_is_fine(&hdr, &buf);
    let _ = assertion_before_a_debug_assert_is_fine(true, 2, &buf);
    let _ = range_pattern_is_a_dispatch_so_is_fine(2, &buf);
    let _ = unsafe { write_through_a_copied_pointer_is_fine(&mut owned, &buf) };
}
