#[macro_export]
macro_rules! vec_of_owned {
    () => (
        vec![]
    );
    ($($x:expr),+ $(,)?) => (
        vec![$($x),+].into_iter().map(|x| x.to_owned()).collect_vec()
    );
}

#[macro_export]
macro_rules! vec_u8 {
    ($($elem:expr),* $(,)?) => {
        vec![$($elem as u8),*]
    };
}
