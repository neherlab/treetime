// An argument named as one parameter must not be bound to another of the same type.

fn resize(width: u32, height: u32) -> u32 {
    width * 2 + height
}

struct SpawnOptions {
    inherit_stdout: bool,
    inherit_stderr: bool,
}

struct Daemon {
    detached: bool,
}

fn spawn(inherit_stdout: bool, inherit_stderr: bool) -> bool {
    inherit_stdout && !inherit_stderr
}

fn place(line: u32, column: u32, offset: u32) -> u32 {
    line + column + offset
}

fn open(path: &str, is_dir: bool, follow: bool) -> usize {
    path.len() + usize::from(is_dir) + usize::from(follow)
}

fn label(name: &str, alias: &str) -> usize {
    name.len() + alias.len()
}

fn trace(index: usize, from_index: usize) -> usize {
    index - from_index
}

fn scale(width: u32, height: f64) -> f64 {
    f64::from(width) * height
}

struct Canvas;

impl Canvas {
    fn blit(&self, src: usize, dst: usize) -> usize {
        dst - src
    }
}

fn under(url: &str, registry: &str) -> bool {
    url.starts_with(registry)
}

fn has_prefix(self_: &[u8], str: &[u8]) -> bool {
    self_.starts_with(str)
}

struct Span(u32, u32);

impl Span {
    fn within(&self, other: &Span) -> bool {
        other.0 <= self.0 && self.1 <= other.1
    }

    fn encloses(&self, inner: &Span) -> bool {
        // Fine: `self` in an argument slot names a position, not a role.
        inner.within(self)
    }
}

fn swapped_pair_is_flagged(width: u32, height: u32) -> u32 {
    // Flagged: both names cross, reported once for the call.
    resize(height, width)
}

fn swapped_fields_are_flagged(opts: &SpawnOptions) -> bool {
    // Flagged: the field names cross the parameter names.
    spawn(opts.inherit_stderr, opts.inherit_stdout)
}

fn one_misbound_is_flagged(column: u32, extra: u32) -> u32 {
    // Flagged: `column` lands in `line` while a `column` parameter exists.
    place(column, extra, 0)
}

fn method_args_are_flagged(c: &Canvas, src: usize, dst: usize) -> usize {
    // Flagged: receiver aside, the two indices are transposed.
    c.blit(dst, src)
}

fn lone_reversal_is_flagged(url: &str, registry: &str) -> bool {
    // Flagged: nothing nearby applies `under` the right way round.
    under(registry, url)
}

fn literal_partner_is_flagged(opts: &SpawnOptions, height: u32) -> u32 {
    // Flagged: a literal in the namesake slot is still the other half of a
    // transposition; `spawn(opts.inherit_stdout, true)` compiles just as well.
    let _ = spawn(true, opts.inherit_stdout);
    resize(height, 0)
}

fn correct_order_is_fine(width: u32, height: u32, opts: &SpawnOptions) -> u32 {
    // Fine: every name sits in its own slot.
    let _ = spawn(opts.inherit_stdout, opts.inherit_stderr);
    resize(width, height)
}

fn prefixed_names_are_fine(dir: bool, is_follow: bool) -> usize {
    // Fine: `dir` is `is_dir` and `is_follow` is `follow` once prefixes go.
    open("p", dir, is_follow)
}

fn different_types_are_fine(height: u32, width: f64) -> f64 {
    // Fine: `scale`'s own `height` is `f64`, so a `u32` named `height` in
    // the `width` slot cannot be the transposed one.
    scale(height, width)
}

fn same_value_twice_is_fine(name: &str) -> usize {
    // Fine: `name` also fills `name`; nothing is transposed.
    label(name, name)
}

fn qualified_param_is_fine(index: usize, next: usize) -> usize {
    // Fine: `from_index` receiving `index` is the recursion's parent, not a swap.
    trace(next, index)
}

fn unnamed_args_are_fine(d: &Daemon) -> bool {
    // Fine: a literal and an unrelated field carry no crossing name.
    spawn(true, d.detached)
}

fn symmetric_pair_is_fine(url: &str, scope_registry: &str) -> bool {
    // Fine: both orders in one condition is an equality test, not a slip.
    !(under(url, scope_registry) && under(scope_registry, url))
}

fn pseudo_receiver_is_fine(str: &[u8]) -> bool {
    // Fine: `self_` is a receiver slot; whatever fills it is the subject.
    has_prefix(str, b"./")
}

fn closures_are_fine(width: u32, height: u32) -> u32 {
    // Fine: a closure's parameters are not a signature anyone reads by name.
    let f = |width: u32, height: u32| width + height;
    f(height, width)
}

fn main() {
    let opts = SpawnOptions {
        inherit_stdout: true,
        inherit_stderr: false,
    };
    let _ = swapped_pair_is_flagged(1, 2);
    let _ = swapped_fields_are_flagged(&opts);
    let _ = one_misbound_is_flagged(1, 2);
    let _ = method_args_are_flagged(&Canvas, 1, 2);
    let _ = correct_order_is_fine(1, 2, &opts);
    let (dir, follow) = (true, opts.inherit_stderr);
    let _ = prefixed_names_are_fine(dir, follow);
    let _ = different_types_are_fine(2, 1.0);
    let _ = same_value_twice_is_fine("n");
    let _ = qualified_param_is_fine(2, 3);
    let _ = unnamed_args_are_fine(&Daemon { detached: false });
    let _ = closures_are_fine(1, 2);
    let _ = lone_reversal_is_flagged("u", "r");
    let _ = literal_partner_is_flagged(&opts, 2);
    let _ = symmetric_pair_is_fine("u", "r");
    let _ = pseudo_receiver_is_fine(b"s");
    let _ = Span(1, 2).encloses(&Span(0, 3));
}
