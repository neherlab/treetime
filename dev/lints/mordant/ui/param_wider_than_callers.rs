// A panicking arm for a variant that no existing call site passes: the
// parameter type is wider than the function's real domain.

enum Shape {
    Circle(u32),
    Square(u32),
    Line,
}

// Flagged: both call sites pass Circle or Square; the Line arm's panic
// cannot fire today and should be a compile error tomorrow.
fn area(s: Shape) -> u32 {
    match s {
        Shape::Circle(r) => 3 * r * r,
        Shape::Square(w) => w * w,
        Shape::Line => unreachable!("lines have no area"),
    }
}

// Fine: a call site passes Line, so the panic arm is genuinely reachable.
fn perimeter(s: Shape) -> u32 {
    match s {
        Shape::Circle(r) => 6 * r,
        Shape::Square(w) => 4 * w,
        Shape::Line => panic!("lines have no perimeter"),
    }
}

// Fine: the argument is not a constructor literal, so the call-site set is
// unknowable and the lint stays silent.
fn diameter(s: Shape) -> u32 {
    match s {
        Shape::Circle(r) => 2 * r,
        _ => unreachable!("only circles have diameters"),
    }
}

fn pick(n: u32) -> Shape {
    if n > 2 { Shape::Circle(n) } else { Shape::Square(n) }
}

impl Shape {
    // Fine: `pick(1).sides()` sends whatever `pick` returned to `self`, so
    // the one path-style call passing Circle proves nothing about Line.
    fn sides(self) -> u32 {
        match self {
            Shape::Circle(_) => 0,
            Shape::Square(_) => 4,
            Shape::Line => unreachable!("lines have no sides"),
        }
    }

    // Flagged: every receiver is a literal, and none of them is Line.
    fn corners(self) -> u32 {
        match self {
            Shape::Circle(_) => 0,
            Shape::Square(_) => 4,
            Shape::Line => unreachable!("lines have no corners"),
        }
    }
}

// Fine: `Shape::Circle` in the first argument position is a constructor
// passed as a function, not a value of `Shape`, so it names no variant; the
// second argument is what `build` matches on, and a call passes `Line`.
fn build(make: fn(u32) -> Shape, seed: Shape) -> Shape {
    match seed {
        Shape::Circle(n) | Shape::Square(n) => make(n),
        Shape::Line => panic!("nothing to build from"),
    }
}

fn main() {
    let _ = area(Shape::Circle(2));
    let _ = area(Shape::Square(3));
    let _ = perimeter(Shape::Circle(2));
    let _ = perimeter(Shape::Line);
    let _ = diameter(pick(4));
    let _ = Shape::sides(Shape::Circle(1));
    let _ = pick(1).sides();
    let _ = Shape::Circle(1).corners();
    let _ = Shape::Square(2).corners();
    let _ = build(Shape::Circle, Shape::Square(1));
    let _ = build(Shape::Square, Shape::Line);
}
