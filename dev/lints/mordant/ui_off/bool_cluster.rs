// The same shape ui/bool_cluster.rs flags; without `bool-cluster-enabled`
// nothing is reported.

struct Refusals {
    too_big: bool,
    exhausted: bool,
    fragmented: bool,
}

fn main() {
    let r = Refusals {
        too_big: false,
        exhausted: false,
        fragmented: true,
    };
    println!("{} {} {}", r.too_big, r.exhausted, r.fragmented);
}
