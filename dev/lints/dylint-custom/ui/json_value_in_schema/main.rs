mod schemars {
    pub trait JsonSchema {}
}

mod serde_json {
    pub enum Value {
        Null,
    }

    pub struct Map<K, V> {
        pub entries: Vec<(K, V)>,
    }
}

use schemars::JsonSchema;
use serde_json::{Map, Value};

pub struct Open {
    pub name: String,
    pub payload: Value,
    pub settings: Map<String, Value>,
    pub items: Vec<Value>,
}

impl JsonSchema for Open {}

pub enum Event {
    Data { payload: Option<Value> },
    Empty,
}

impl JsonSchema for Event {}

pub struct Typed {
    pub name: String,
}

impl JsonSchema for Typed {}

pub struct NotInSchema {
    pub payload: Value,
}

pub struct JsonValue(pub Value);

impl JsonSchema for JsonValue {}

fn main() {
    let _ = (
        Open { name: String::new(), payload: Value::Null, settings: Map { entries: vec![] }, items: vec![] },
        Event::Data { payload: None },
        Event::Empty,
        Typed { name: String::new() },
        NotInSchema { payload: Value::Null },
        JsonValue(Value::Null),
    );
}
