// The shape ui/unchecked_input_len.rs flags first; without
// `unchecked-input-len-enabled` nothing is reported.

fn use_before_check(len: usize, buf: &[u8]) -> Option<&[u8]> {
    let (head, _) = buf.split_at(len);
    if len > buf.len() {
        return None;
    }
    Some(head)
}

fn main() {
    let buf = [0u8; 8];
    let _ = use_before_check(2, &buf);
}
