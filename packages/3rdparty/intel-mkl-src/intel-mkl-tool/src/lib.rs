//! Library files in Intel MKL 2020.1 for Linux
//! --------------------------------------------
//!
//! ### MKL Core
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_core                |                  libmkl_core.so|                  libmkl_core.a|
//! |mkl_def                 |                   libmkl_def.so|                             - |
//! |mkl_rt                  |                    libmkl_rt.so|                             - |
//! |mkl_avx                 |                   libmkl_avx.so|                             - |
//! |mkl_avx2                |                  libmkl_avx2.so|                             - |
//! |mkl_avx512              |                libmkl_avx512.so|                             - |
//! |mkl_avx512_mic          |            libmkl_avx512_mic.so|                             - |
//! |mkl_mc                  |                    libmkl_mc.so|                             - |
//! |mkl_mc3                 |                   libmkl_mc3.so|                             - |
//!
//! ### Vector Math library
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_vml_def             |               libmkl_vml_def.so|                             - |
//! |mkl_vml_avx             |               libmkl_vml_avx.so|                             - |
//! |mkl_vml_avx2            |              libmkl_vml_avx2.so|                             - |
//! |mkl_vml_avx512          |            libmkl_vml_avx512.so|                             - |
//! |mkl_vml_avx512_mic      |        libmkl_vml_avx512_mic.so|                             - |
//! |mkl_vml_mc              |                libmkl_vml_mc.so|                             - |
//! |mkl_vml_mc2             |               libmkl_vml_mc2.so|                             - |
//! |mkl_vml_mc3             |               libmkl_vml_mc3.so|                             - |
//! |mkl_vml_cmpt            |              libmkl_vml_cmpt.so|                             - |
//!
//! ### Intel OpenMP
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |iomp5                   |                     libiomp5.so|                     libiomp5.a|
//! |iomp5_db                |                  libiomp5_db.so|                             - |
//! |iompstubs5              |                libiompstubs5.so|                libiompstubs5.a|
//!
//! ### Threading switch
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_sequential          |            libmkl_sequential.so|            libmkl_sequential.a|
//! |mkl_gnu_thread          |            libmkl_gnu_thread.so|            libmkl_gnu_thread.a|
//! |mkl_intel_thread        |          libmkl_intel_thread.so|          libmkl_intel_thread.a|
//! |mkl_pgi_thread          |            libmkl_pgi_thread.so|            libmkl_pgi_thread.a|
//! |mkl_tbb_thread          |            libmkl_tbb_thread.so|            libmkl_tbb_thread.a|
//!
//! ### LP64/ILP64 switch for Intel and GCC Fortran compilers
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_intel_ilp64         |           libmkl_intel_ilp64.so|           libmkl_intel_ilp64.a|
//! |mkl_intel_lp64          |            libmkl_intel_lp64.so|            libmkl_intel_lp64.a|
//! |mkl_gf_ilp64            |              libmkl_gf_ilp64.so|              libmkl_gf_ilp64.a|
//! |mkl_gf_lp64             |               libmkl_gf_lp64.so|               libmkl_gf_lp64.a|
//!
//! ### Fortran 95 interface
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_blas95_ilp64        |                              - |          libmkl_blas95_ilp64.a|
//! |mkl_blas95_lp64         |                              - |           libmkl_blas95_lp64.a|
//! |mkl_lapack95_ilp64      |                              - |        libmkl_lapack95_ilp64.a|
//! |mkl_lapack95_lp64       |                              - |         libmkl_lapack95_lp64.a|
//!
//! ### ScaLAPACK
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_scalapack_ilp64     |       libmkl_scalapack_ilp64.so|       libmkl_scalapack_ilp64.a|
//! |mkl_scalapack_lp64      |        libmkl_scalapack_lp64.so|        libmkl_scalapack_lp64.a|
//!
//! ### BLACS
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_blacs_intelmpi_ilp64|  libmkl_blacs_intelmpi_ilp64.so|  libmkl_blacs_intelmpi_ilp64.a|
//! |mkl_blacs_intelmpi_lp64 |   libmkl_blacs_intelmpi_lp64.so|   libmkl_blacs_intelmpi_lp64.a|
//! |mkl_blacs_openmpi_ilp64 |   libmkl_blacs_openmpi_ilp64.so|   libmkl_blacs_openmpi_ilp64.a|
//! |mkl_blacs_openmpi_lp64  |    libmkl_blacs_openmpi_lp64.so|    libmkl_blacs_openmpi_lp64.a|
//! |mkl_blacs_sgimpt_ilp64  |    libmkl_blacs_sgimpt_ilp64.so|    libmkl_blacs_sgimpt_ilp64.a|
//! |mkl_blacs_sgimpt_lp64   |     libmkl_blacs_sgimpt_lp64.so|     libmkl_blacs_sgimpt_lp64.a|
//!
//! ### FFT
//!
//! |                        |                         shared |                       static  |
//! |:-----------------------|:-------------------------------|:------------------------------|
//! |mkl_cdft_core           |             libmkl_cdft_core.so|             libmkl_cdft_core.a|
//!

#![cfg_attr(not(feature = "archive"), allow(dead_code))]

use anyhow::*;
use std::path::*;

mod config;
mod entry;

#[cfg(feature = "archive")]
mod download;
#[cfg(feature = "archive")]
mod package;

pub use config::*;
pub use entry::*;

const S3_ADDR: &'static str = "https://s3-ap-northeast-1.amazonaws.com/rust-intel-mkl";

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
mod mkl {
    pub const OS: &str = "linux";
    pub const EXTENSION_STATIC: &'static str = "a";
    pub const EXTENSION_SHARED: &'static str = "so";
    pub const PREFIX: &'static str = "lib";
    pub const VERSION_YEAR: u32 = 2020;
    pub const VERSION_UPDATE: u32 = 1;
}

#[cfg(all(target_os = "macos", target_arch = "x86_64"))]
mod mkl {
    pub const OS: &str = "macos";
    pub const EXTENSION_STATIC: &'static str = "a";
    pub const EXTENSION_SHARED: &'static str = "dylib";
    pub const PREFIX: &'static str = "lib";
    pub const VERSION_YEAR: u32 = 2019;
    pub const VERSION_UPDATE: u32 = 3;
}

#[cfg(all(target_os = "windows", target_arch = "x86_64"))]
mod mkl {
    pub const OS: &str = "windows";
    pub const EXTENSION_STATIC: &'static str = "lib";
    pub const EXTENSION_SHARED: &'static str = "lib";
    pub const PREFIX: &'static str = "";
    pub const VERSION_YEAR: u32 = 2020;
    pub const VERSION_UPDATE: u32 = 1;
}

fn s3_addr() -> String {
    format!(
        "{}/{}/{}.{}",
        S3_ADDR,
        mkl::OS,
        mkl::VERSION_YEAR,
        mkl::VERSION_UPDATE
    )
}

pub fn xdg_home_path() -> PathBuf {
    dirs::data_local_dir().unwrap().join(format!(
        "intel-mkl-tool/{}.{}",
        mkl::VERSION_YEAR,
        mkl::VERSION_UPDATE
    ))
}
