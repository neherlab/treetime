fn main() {
    let x: u64 = 1;
    let _ = x * 1000;
    let _ = x * 7;
}

#[test]
fn test_function() {
    let x: u64 = 1;
    let _ = x * 1001;
}

#[cfg(test)]
mod tests {
    fn helper() -> u64 {
        let x: u64 = 1;
        x * 1002
    }

    #[test]
    fn nested_test_function() {
        let _ = helper() * 1003;
    }
}
