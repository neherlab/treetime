// The shape ui/generic_body_not_generic.rs flags first. Without
// `generic-body-not-generic-enabled` nothing is reported.

pub fn checksum<B: AsRef<[u8]>>(bytes: B) -> u32 {
    let src = bytes.as_ref();
    let mut acc = 17u32;
    let mut run = 0u32;
    for byte in src {
        let v = u32::from(*byte);
        acc = acc.wrapping_mul(31).wrapping_add(v);
        if v & 1 == 0 {
            run += 1;
        } else {
            run = 0;
        }
        acc ^= run << 3;
    }
    (acc ^ (src.len() as u32)).rotate_left(7)
}

fn main() {
    let _ = checksum("abc");
    let _ = checksum(vec![1u8, 2, 3]);
    let _ = checksum([9u8; 4]);
}
