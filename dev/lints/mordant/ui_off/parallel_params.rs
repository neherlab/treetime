// The same shape ui/parallel_params.rs flags; without `parallel-params-enabled`
// nothing is reported.

fn decode(src: &[u8], w: u32, h: u32) -> Vec<u8> {
    let out = scale(src, w, h, 2);
    let _ = checksum(w, h, &out);
    out
}

fn scale(src: &[u8], w: u32, h: u32, factor: u32) -> Vec<u8> {
    src.iter()
        .cycle()
        .take((w * h * factor) as usize)
        .copied()
        .collect()
}

fn checksum(w: u32, h: u32, px: &[u8]) -> u64 {
    u64::from(w) * u64::from(h) + px.len() as u64
}

fn main() {
    let _ = decode(&[1, 2, 3], 4, 4);
}
