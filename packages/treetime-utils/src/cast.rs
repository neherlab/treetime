use num_traits::ToPrimitive;

pub const F64_EXACT_INTEGER_LIMIT: u64 = 1 << f64::MANTISSA_DIGITS;

pub trait ToF64Exact: Copy {
  fn to_f64_exact(self) -> f64;
}

macro_rules! impl_to_f64_exact_lossless {
  ($($int:ty),*) => {$(
    impl ToF64Exact for $int {
      #[inline]
      fn to_f64_exact(self) -> f64 {
        f64::from(self)
      }
    }
  )*};
}

macro_rules! impl_to_f64_exact_bounded {
  ($($int:ty => |$value:ident| $magnitude:expr),*) => {$(
    impl ToF64Exact for $int {
      #[inline]
      #[expect(
        clippy::as_conversions,
        reason = "the debug assertion bounds the magnitude by 2^53, below which every integer converts to f64 exactly"
      )]
      fn to_f64_exact(self) -> f64 {
        let $value = self;
        debug_assert!(
          u64::try_from($magnitude).is_ok_and(|magnitude| magnitude <= F64_EXACT_INTEGER_LIMIT),
          "integer {self} has no exact f64 representation",
        );
        self as f64
      }
    }
  )*};
}

impl_to_f64_exact_lossless!(u8, u16, u32, i8, i16, i32);
impl_to_f64_exact_bounded!(
  usize => |value| value,
  u64 => |value| value,
  isize => |value| value.unsigned_abs(),
  i64 => |value| value.unsigned_abs()
);

const _: () = assert!(usize::BITS >= u32::BITS && usize::BITS <= u64::BITS);

#[expect(
  clippy::as_conversions,
  reason = "the compile-time assertion above guarantees that usize is at most 64 bits wide"
)]
pub const fn usize_to_u64(value: usize) -> u64 {
  value as u64
}

#[expect(
  clippy::as_conversions,
  reason = "the compile-time assertion above guarantees that usize is at least 32 bits wide"
)]
pub const fn u32_to_usize(value: u32) -> usize {
  value as usize
}

pub fn f64_to_usize_saturating(value: f64) -> usize {
  value.to_usize().unwrap_or(if value > 0.0 { usize::MAX } else { 0 })
}

pub fn f64_to_i32_saturating(value: f64) -> i32 {
  value.to_i32().unwrap_or(if value > 0.0 {
    i32::MAX
  } else if value < 0.0 {
    i32::MIN
  } else {
    0
  })
}
