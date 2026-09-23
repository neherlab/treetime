fn main() {
    let _ = String::from(Adjacent::A);
    let _ = String::from(Separated::A);
}

pub enum Adjacent {
    A,
}

impl From<Adjacent> for String {
    fn from(value: Adjacent) -> Self {
        match value {
            Adjacent::A => "adjacent".to_owned(),
        }
    }
}

pub enum Separated {
    A,
}

pub struct Unrelated;

impl From<Separated> for String {
    fn from(value: Separated) -> Self {
        match value {
            Separated::A => "separated".to_owned(),
        }
    }
}
