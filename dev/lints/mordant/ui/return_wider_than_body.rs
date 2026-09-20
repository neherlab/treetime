// A panicking arm for a variant the callee provably never constructs: the
// return type promises more than the function delivers.

enum Token {
    Word(u32),
    Space,
    Eof,
}

// Every return is a constructor literal, and none of them is Eof.
fn next_token(n: u32) -> Token {
    if n == 0 {
        return Token::Space;
    }
    match n {
        1 => Token::Word(1),
        _ => Token::Word(n),
    }
}

// Flagged: next_token never constructs Eof.
fn count(n: u32) -> u32 {
    match next_token(n) {
        Token::Word(w) => w,
        Token::Space => 0,
        Token::Eof => unreachable!("the tokenizer never yields Eof here"),
    }
}

// Flagged too: the MIR dataflow follows the local binding and both branch
// assignments, so indirection does not make the set unknowable.
fn via_local(n: u32) -> Token {
    let t = if n > 0 { Token::Word(n) } else { Token::Space };
    t
}

fn count_via_local(n: u32) -> u32 {
    match via_local(n) {
        Token::Word(w) => w,
        Token::Space => 0,
        Token::Eof => unreachable!("never yielded by via_local"),
    }
}

// Fine: this producer can return Eof, so the panic arm is reachable.
fn next_or_eof(n: u32) -> Token {
    if n > 100 { Token::Eof } else { Token::Word(n) }
}

fn count_or_eof(n: u32) -> u32 {
    match next_or_eof(n) {
        Token::Word(w) => w,
        Token::Space => 0,
        Token::Eof => panic!("stream ended"),
    }
}

// Fine: one return position is not a constructor literal, so the set is
// unknowable and the lint stays silent about this producer.
fn passthrough(t: Token) -> Token {
    match t {
        Token::Word(w) => Token::Word(w + 1),
        other => other,
    }
}

fn count_passthrough(t: Token) -> u32 {
    match passthrough(t) {
        Token::Word(w) => w,
        Token::Space => 0,
        Token::Eof => panic!("no eof expected"),
    }
}

fn main() {
    let _ = count(1) + count_or_eof(2) + count_passthrough(Token::Space) + count_via_local(3);
}
