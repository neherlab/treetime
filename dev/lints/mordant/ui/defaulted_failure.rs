// A callee that rejects some of its argument, and a caller that replaces the
// rejection with a fixed value and carries on. The `ui` config lists
// `from_str_radix` and `listed_by_config` under `defaulted-failure-callees`
// and `Pending` under `defaulted-failure-ignored-errors`.

enum ParseError {
    Empty,
    TooBig,
}

fn describe(e: &ParseError) -> &'static str {
    match e {
        ParseError::Empty => "empty",
        ParseError::TooBig => "too big",
    }
}

// Rejects on the argument: a validator.
fn parse_port(s: &str) -> Result<u16, ParseError> {
    if s.is_empty() {
        return Err(ParseError::Empty);
    }
    let n: u32 = s.len() as u32 * 1000;
    if n > u16::MAX as u32 {
        return Err(ParseError::TooBig);
    }
    Ok(n as u16)
}

// Rejects on the argument, but as an `Option`: absent, not refused.
fn parse_flag(s: &str) -> Option<bool> {
    match s {
        "on" => Some(true),
        "off" => Some(false),
        _ => None,
    }
}

// Rejects through `?` on a check of the argument.
fn parse_pair(s: &str) -> Result<(u16, u16), ParseError> {
    let (a, b) = s.split_once(':').ok_or(ParseError::Empty)?;
    Ok((parse_port(a)?, parse_port(b)?))
}

struct Header {
    limit: u32,
    pos: usize,
}

impl Header {
    // A method whose failure the argument decides (together with `self`).
    fn parse_limit(&self, s: &str) -> Result<u32, ParseError> {
        let n = s.len() as u32;
        if n > self.limit {
            return Err(ParseError::TooBig);
        }
        Ok(n)
    }

    // A method whose failure only its own state decides.
    fn next_field(&mut self, s: &str) -> Result<u32, ParseError> {
        if self.pos >= 4 {
            return Err(ParseError::Empty);
        }
        self.pos += 1;
        Ok(s.len() as u32)
    }

    // A `bool` answer.
    fn is_open(&self, s: &str) -> Result<bool, ParseError> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }
        Ok(s.len() < 3)
    }

    // The argument decides an early success; only its own state decides the
    // failure behind that.
    fn send(&self, s: &str) -> Result<u32, ParseError> {
        if s.is_empty() {
            return Ok(0);
        }
        if self.pos >= 4 {
            return Err(ParseError::Empty);
        }
        Ok(s.len() as u32)
    }

    // The other way round: the argument decides the failure, once its own
    // state has let it get that far.
    fn send_checked(&self, s: &str) -> Result<u32, ParseError> {
        if self.pos >= 4 {
            if s.is_empty() {
                return Err(ParseError::Empty);
            }
        }
        Ok(s.len() as u32)
    }
}

// Rejects on the argument, with nothing to hand back on success.
fn check_port(s: &str) -> Result<(), ParseError> {
    parse_port(s)?;
    Ok(())
}

// An error that is recorded elsewhere by the time it is returned; the config
// lists it, so defaulting it hides nothing.
struct Pending;

fn schedule(s: &str, log: &mut Vec<String>) -> Result<u32, Pending> {
    if s.is_empty() {
        log.push("empty".to_string());
        return Err(Pending);
    }
    Ok(s.len() as u32)
}

// Fails only because a resource refused: never a verdict on the path.
fn read_size(path: &str) -> Result<u64, std::io::Error> {
    std::fs::metadata(path).map(|m| m.len())
}

// Fails on state it was not handed, not on the argument.
fn next_slot(tag: &str) -> Option<u32> {
    let free = std::env::var_os("SLOTS").is_some();
    if !free {
        return None;
    }
    Some(tag.len() as u32)
}

// The same, as a `Result`, behind an early success the argument decides.
fn slot_for(tag: &str) -> Result<u32, ParseError> {
    if tag.is_empty() {
        return Ok(0);
    }
    if std::env::var_os("SLOTS").is_none() {
        return Err(ParseError::Empty);
    }
    Ok(tag.len() as u32)
}

// Cannot fail.
fn always(s: &str) -> Option<usize> {
    Some(s.len())
}

// Its failure is built by a combinator, so the body shows no check of its
// own; only the config makes it visible.
fn listed_by_config(s: &str) -> Option<u32> {
    s.parse().ok()
}

fn unwrap_or_literal_is_flagged(s: &str) -> u16 {
    parse_port(s).unwrap_or(0)
}

fn unwrap_or_default_is_flagged(s: &str) -> u16 {
    parse_port(s).unwrap_or_default()
}

fn unwrap_or_else_const_is_flagged(s: &str) -> u16 {
    parse_port(s).unwrap_or_else(|_| u16::MAX)
}

fn unwrap_or_else_default_call_is_flagged(s: &str) -> u16 {
    parse_port(s).unwrap_or_else(|_| Default::default())
}

fn ok_then_unwrap_or_is_flagged(s: &str) -> u16 {
    parse_port(s).ok().unwrap_or(80)
}

fn bool_defaulted_to_true_is_flagged(h: &Header, s: &str) -> bool {
    h.is_open(s).unwrap_or(true)
}

fn question_mark_callee_is_flagged(s: &str) -> (u16, u16) {
    parse_pair(s).unwrap_or((0, 0))
}

fn method_callee_is_flagged(h: &Header, s: &str) -> u32 {
    h.parse_limit(s).unwrap_or(0)
}

fn let_else_ok_unit_is_flagged(s: &str, out: &mut Vec<u16>) -> Result<(), ParseError> {
    let Ok(port) = parse_port(s) else {
        return Ok(());
    };
    out.push(port);
    Ok(())
}

fn let_else_bare_return_is_flagged(s: &str, out: &mut Vec<u16>) {
    let Ok(port) = parse_port(s) else { return };
    out.push(port);
}

fn let_else_some_of_ok_unit_is_flagged(s: &str, out: &mut Vec<u16>) -> Result<(), ParseError> {
    let Some(port) = parse_port(s).ok() else {
        return Ok(());
    };
    out.push(port);
    Ok(())
}

fn let_else_some_of_ok_bare_return_is_flagged(s: &str, out: &mut Vec<u16>) {
    let Some(port) = parse_port(s).ok() else { return };
    out.push(port);
}

fn argument_check_behind_receiver_check_is_flagged(h: &Header, s: &str) -> u32 {
    h.send_checked(s).unwrap_or(0)
}

fn unit_default_in_value_position_is_flagged(s: &str) {
    check_port(s).unwrap_or_default()
}

fn listed_callee_is_flagged(s: &str) -> u32 {
    listed_by_config(s).unwrap_or(0)
}

fn listed_foreign_callee_is_flagged(s: &str) -> u32 {
    u32::from_str_radix(s, 16).unwrap_or(0)
}

fn option_callee_is_fine(s: &str) -> bool {
    parse_flag(s).unwrap_or(true)
}

fn option_callee_default_fn_is_fine(s: &str) -> bool {
    parse_flag(s).unwrap_or_else(Default::default)
}

fn receiver_decided_is_fine(h: &mut Header, s: &str) -> u32 {
    h.next_field(s).unwrap_or(0)
}

fn bool_defaulted_to_false_is_fine(h: &Header, s: &str) -> bool {
    h.is_open(s).unwrap_or(false)
}

fn bool_defaulted_by_default_is_fine(h: &Header, s: &str) -> bool {
    h.is_open(s).unwrap_or_default()
}

fn bool_defaulted_to_false_by_thunk_is_fine(h: &Header, s: &str) -> bool {
    h.is_open(s).unwrap_or_else(|_| false)
}

fn bool_defaulted_by_default_call_in_thunk_is_fine(h: &Header, s: &str) -> bool {
    h.is_open(s).unwrap_or_else(|_| Default::default())
}

fn bool_defaulted_by_default_fn_is_fine(h: &Header, s: &str) -> bool {
    h.is_open(s).ok().unwrap_or_else(Default::default)
}

fn early_success_before_receiver_check_is_fine(h: &Header, s: &str) -> u32 {
    h.send(s).unwrap_or(0)
}

fn early_success_before_state_check_is_fine(tag: &str) -> u32 {
    slot_for(tag).unwrap_or(0)
}

fn statement_position_default_is_a_discard(s: &str) {
    check_port(s).unwrap_or_default();
}

fn ignored_error_is_fine(s: &str, log: &mut Vec<String>) -> u32 {
    schedule(s, log).unwrap_or(0)
}

fn foreign_callee_is_fine(s: &str) -> u32 {
    s.parse().unwrap_or(0)
}

fn resource_failure_is_fine(path: &str) -> u64 {
    read_size(path).unwrap_or(0)
}

fn state_failure_is_fine(tag: &str) -> u32 {
    next_slot(tag).unwrap_or(0)
}

fn infallible_is_fine(s: &str) -> usize {
    always(s).unwrap_or(0)
}

fn computed_fallback_is_fine(s: &str, prev: u16) -> u16 {
    parse_port(s).unwrap_or(prev)
}

fn computed_thunk_is_fine(s: &str, prev: u16) -> u16 {
    parse_port(s).unwrap_or_else(|_| prev + 1)
}

fn handled_is_fine(s: &str) -> u16 {
    match parse_port(s) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("{}", describe(&e));
            0
        }
    }
}

fn propagated_is_fine(s: &str) -> Result<u16, ParseError> {
    let p = parse_port(s)?;
    Ok(p)
}

#[allow(unknown_lints, discarded_error)]
fn statement_ok_is_another_lints(s: &str) {
    parse_port(s).ok();
}

fn let_else_some_is_fine(s: &str, out: &mut Vec<bool>) {
    let Some(flag) = parse_flag(s) else { return };
    out.push(flag);
}

fn let_else_reporting_failure_is_fine(s: &str, out: &mut Vec<u16>) -> bool {
    let Ok(port) = parse_port(s) else {
        return false;
    };
    out.push(port);
    true
}

fn let_else_returning_none_is_fine(s: &str) -> Option<u16> {
    let Ok(port) = parse_port(s) else {
        return None;
    };
    Some(port)
}

fn let_else_doing_work_is_fine(s: &str, log: &mut Vec<String>) {
    let Ok(port) = parse_port(s) else {
        log.push(format!("bad port {s}"));
        return;
    };
    log.push(port.to_string());
}

fn main() {
    let h = Header { limit: 4, pos: 0 };
    let mut hm = Header { limit: 4, pos: 0 };
    let mut ports = Vec::new();
    let mut flags = Vec::new();
    let mut log = Vec::new();
    let _ = unwrap_or_literal_is_flagged("a");
    let _ = unwrap_or_default_is_flagged("a");
    let _ = unwrap_or_else_const_is_flagged("a");
    let _ = unwrap_or_else_default_call_is_flagged("a");
    let _ = ok_then_unwrap_or_is_flagged("a");
    let _ = bool_defaulted_to_true_is_flagged(&h, "a");
    let _ = question_mark_callee_is_flagged("a:b");
    let _ = method_callee_is_flagged(&h, "a");
    let _ = let_else_ok_unit_is_flagged("a", &mut ports);
    let_else_bare_return_is_flagged("a", &mut ports);
    let _ = let_else_some_of_ok_unit_is_flagged("a", &mut ports);
    let_else_some_of_ok_bare_return_is_flagged("a", &mut ports);
    let _ = argument_check_behind_receiver_check_is_flagged(&h, "a");
    unit_default_in_value_position_is_flagged("a");
    let _ = listed_callee_is_flagged("1");
    let _ = listed_foreign_callee_is_flagged("ff");
    let _ = option_callee_is_fine("on");
    let _ = option_callee_default_fn_is_fine("on");
    let _ = receiver_decided_is_fine(&mut hm, "a");
    let _ = bool_defaulted_to_false_is_fine(&h, "a");
    let _ = bool_defaulted_by_default_is_fine(&h, "a");
    let _ = bool_defaulted_to_false_by_thunk_is_fine(&h, "a");
    let _ = bool_defaulted_by_default_call_in_thunk_is_fine(&h, "a");
    let _ = bool_defaulted_by_default_fn_is_fine(&h, "a");
    let _ = early_success_before_receiver_check_is_fine(&h, "a");
    let _ = early_success_before_state_check_is_fine("a");
    statement_position_default_is_a_discard("a");
    let _ = ignored_error_is_fine("a", &mut log);
    let _ = foreign_callee_is_fine("1");
    let _ = resource_failure_is_fine("/");
    let _ = state_failure_is_fine("a");
    let _ = infallible_is_fine("a");
    let _ = computed_fallback_is_fine("a", 1);
    let _ = computed_thunk_is_fine("a", 1);
    let _ = handled_is_fine("a");
    let _ = propagated_is_fine("a");
    statement_ok_is_another_lints("a");
    let_else_some_is_fine("on", &mut flags);
    let _ = let_else_reporting_failure_is_fine("a", &mut ports);
    let _ = let_else_returning_none_is_fine("a");
    let_else_doing_work_is_fine("a", &mut log);
}
