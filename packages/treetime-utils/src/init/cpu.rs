use std::io::{self, Write};

pub const CPU_FEATURES_MESSAGE: &str = cpu_features_message();

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    debug_remnants,
    reason = "the message reaches the user before the logger exists, because the check runs first in the constructor of the binary"
  )
)]
pub fn require_cpu_features() {
  if !cpu_has_build_features() {
    let code = if writeln!(io::stderr().lock(), "{CPU_FEATURES_MESSAGE}").is_ok() {
      1
    } else {
      2
    };
    exit_now(code);
  }
}

#[cfg(unix)]
#[expect(
  unsafe_code,
  reason = "the check runs in a constructor before main, where exit() would also run the destructors of constructors that have not run yet; _exit ends the process without them"
)]
fn exit_now(code: i32) -> ! {
  // SAFETY: _exit has no preconditions and does not return.
  unsafe { libc::_exit(code) }
}

#[cfg(not(unix))]
#[expect(
  clippy::exit,
  reason = "the program stops before code of the build CPU level runs, so no error can propagate to main"
)]
fn exit_now(code: i32) -> ! {
  std::process::exit(code)
}

const fn cpu_features_message() -> &'static str {
  if cfg!(all(target_arch = "x86_64", target_feature = "avx2")) {
    "TreeTime needs a 64-bit x86 CPU with AVX2, FMA, and BMI2, such as Intel Haswell (2013), AMD Excavator (2015), or newer. This CPU lacks some of these features."
  } else if cfg!(all(target_arch = "x86_64", target_feature = "sse4.2")) {
    "TreeTime needs a 64-bit x86 CPU with SSE4.2, POPCNT, and CMPXCHG16B, such as Intel Nehalem (2008), AMD Bulldozer (2011), or newer. This CPU lacks some of these features."
  } else if cfg!(target_arch = "aarch64") {
    "TreeTime needs a 64-bit Arm CPU with the Armv8.2 features LSE atomics, RDM, and CRC32, such as Arm Cortex-A55, Cortex-A76, Neoverse N1, or newer. This CPU lacks some of these features."
  } else {
    "TreeTime does not support this CPU."
  }
}

#[cfg(target_arch = "x86_64")]
const X86_REQUIRED: [u32; 4] = [
  bit(cfg!(target_feature = "sse3"), 0)
    | bit(cfg!(target_feature = "pclmulqdq"), 1)
    | bit(cfg!(target_feature = "ssse3"), 9)
    | bit(cfg!(target_feature = "fma"), 12)
    | bit(cfg!(target_feature = "cmpxchg16b"), 13)
    | bit(cfg!(target_feature = "sse4.1"), 19)
    | bit(cfg!(target_feature = "sse4.2"), 20)
    | bit(cfg!(target_feature = "movbe"), 22)
    | bit(cfg!(target_feature = "popcnt"), 23)
    | bit(cfg!(target_feature = "xsave"), 26)
    | bit(cfg!(target_feature = "avx"), 27)
    | bit(cfg!(target_feature = "avx"), 28)
    | bit(cfg!(target_feature = "f16c"), 29)
    | bit(cfg!(target_feature = "rdrand"), 30),
  bit(cfg!(target_feature = "bmi1"), 3) | bit(cfg!(target_feature = "avx2"), 5) | bit(cfg!(target_feature = "bmi2"), 8),
  bit(cfg!(target_feature = "lahfsahf"), 0) | bit(cfg!(target_feature = "lzcnt"), 5),
  bit(cfg!(target_feature = "avx"), 1) | bit(cfg!(target_feature = "avx"), 2),
];

#[cfg(target_arch = "x86_64")]
pub fn cpu_has_build_features() -> bool {
  x86::features_present(&X86_REQUIRED) != 0
}

#[cfg(target_arch = "x86_64")]
mod x86 {
  #![cfg_attr(
    not(dylint_lib = "custom"),
    expect(
      unsafe_code,
      reason = "the naked function below is unsafe assembly; the nightly toolchain of the lint libraries does not report the unsafe naked attribute"
    )
  )]
  #![expect(
    clippy::inline_asm_x86_intel_syntax,
    reason = "assembly, so that the compiler cannot emit an instruction of the build CPU level before the check has run"
  )]
  #![cfg_attr(
    not(dylint_lib = "custom"),
    expect(
      unsafe_code,
      reason = "the function follows the sysv64 convention, reads the four words behind its argument, saves rbx, which cpuid overwrites, and returns 1 when every required bit of CPUID leaf 1 ECX, leaf 7 EBX, leaf 0x80000001 ECX, and XCR0 is set; the older nightly of the dylint driver does not lint the unsafe naked attribute"
    )
  )]

  #[unsafe(naked)]
  pub(super) extern "sysv64" fn features_present(required: &[u32; 4]) -> u32 {
    core::arch::naked_asm!(
      "push rbx",
      "mov r8, rdi",
      "xor eax, eax",
      "cpuid",
      "mov r9d, eax",
      "mov eax, 1",
      "cpuid",
      "mov eax, ecx",
      "and eax, dword ptr [r8]",
      "cmp eax, dword ptr [r8]",
      "jne 4f",
      "mov eax, dword ptr [r8 + 12]",
      "test eax, eax",
      "jz 2f",
      "xor ecx, ecx",
      "xgetbv",
      "and eax, dword ptr [r8 + 12]",
      "cmp eax, dword ptr [r8 + 12]",
      "jne 4f",
      "2:",
      "mov eax, dword ptr [r8 + 4]",
      "test eax, eax",
      "jz 3f",
      "cmp r9d, 7",
      "jb 4f",
      "mov eax, 7",
      "xor ecx, ecx",
      "cpuid",
      "mov eax, ebx",
      "and eax, dword ptr [r8 + 4]",
      "cmp eax, dword ptr [r8 + 4]",
      "jne 4f",
      "3:",
      "mov eax, 0x80000001",
      "cpuid",
      "mov eax, ecx",
      "and eax, dword ptr [r8 + 8]",
      "cmp eax, dword ptr [r8 + 8]",
      "jne 4f",
      "mov eax, 1",
      "pop rbx",
      "ret",
      "4:",
      "xor eax, eax",
      "pop rbx",
      "ret",
    )
  }
}

#[cfg(all(target_arch = "aarch64", target_os = "linux"))]
#[expect(
  unsafe_code,
  reason = "getauxval is the only source of the hardware capabilities of the kernel"
)]
pub fn cpu_has_build_features() -> bool {
  let required = u64::from(
    bit(cfg!(target_feature = "crc"), 7)
      | bit(cfg!(target_feature = "lse"), 8)
      | bit(cfg!(target_feature = "rdm"), 12)
      | bit(cfg!(target_feature = "dpb"), 16),
  );
  // SAFETY: getauxval has no preconditions; AT_HWCAP is always present on Linux.
  let hwcap = unsafe { libc::getauxval(libc::AT_HWCAP) };
  hwcap & required == required
}

#[cfg(not(any(target_arch = "x86_64", all(target_arch = "aarch64", target_os = "linux"))))]
pub fn cpu_has_build_features() -> bool {
  true
}

#[cfg(any(target_arch = "x86_64", all(target_arch = "aarch64", target_os = "linux")))]
const fn bit(enabled: bool, position: u32) -> u32 {
  if enabled { 1 << position } else { 0 }
}
