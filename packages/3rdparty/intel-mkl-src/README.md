# intel-mkl-src

|crate         | crate.io                                                                                           | description                                                           |
|:-------------|:---------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------|
|intel-mkl-src | [![Crate](http://meritbadge.herokuapp.com/intel-mkl-src)](https://crates.io/crates/intel-mkl-src)  | Source crate for Intel-MKL                                            |
|intel-mkl-sys | [![Crate](http://meritbadge.herokuapp.com/intel-mkl-sys)](https://crates.io/crates/intel-mkl-sys)  | FFI for Intel-MKL [vector math][VM], and [statistical functions][VSL] |
|intel-mkl-tool| [![Crate](http://meritbadge.herokuapp.com/intel-mkl-tool)](https://crates.io/crates/intel-mkl-tool)| CLI utility for redistributing Intel-MKL                              |

Redistribution of Intel MKL as a crate. Tested on Linux, macOS, and Windows (since 0.4.0)

[VM]:  https://software.intel.com/en-us/mkl-developer-reference-c-vector-mathematical-functions
[VSL]: https://software.intel.com/en-us/mkl-developer-reference-c-statistical-functions

## Supported features

- `mkl-*-*-*` features specify which MKL to be linked
  - `static` means MKL will be linked statically, and `dynamic` means MKL will be linked dynamically
  - `lp64` means 32-bit integer interface, `ilp64` means 64-bit integer interface
  - `iomp` means MKL uses Intel OpenMP, `seq` means sequential execution, e.g. no parallelization
    - OpenMP is not supported for Windows currently [#46](https://github.com/rust-math/intel-mkl-src/issues/46)
  - default is `mkl-static-ilp64-seq`, and you must choose one of them.
  - macOS is not supported [#42](https://github.com/rust-math/intel-mkl-src/issues/42)

    | feature name           | Linux              | macOS              | Windows            |
    |:-----------------------|:------------------:|:------------------:|:------------------:|
    | mkl-static-lp64-iomp   | :heavy_check_mark: | -                  | -                  |
    | mkl-static-lp64-seq    | :heavy_check_mark: | -                  | :heavy_check_mark: |
    | mkl-static-ilp64-iomp  | :heavy_check_mark: | -                  | -                  |
    | mkl-static-ilp64-seq   | :heavy_check_mark: | -                  | :heavy_check_mark: |
    | mkl-dynamic-lp64-iomp  | :heavy_check_mark: | :heavy_check_mark: | -                  |
    | mkl-dynamic-lp64-seq   | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
    | mkl-dynamic-ilp64-iomp | :heavy_check_mark: | :heavy_check_mark: | -                  |
    | mkl-dynamic-ilp64-seq  | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |

- `download` feature enables downloading MKL archive managed by this project from AWS S3 (default ON)

## Usage

This crate is a `*-src` crate. This downloads and link Intel MKL, but does not introduce any symbols.
Please use `blas-sys`, `lapack-sys`, or `fftw-sys` to use BLAS, LAPACK, FFTW interface of MKL, e.g.

```toml
[dependencies]
fftw-sys = { version = "0.4", features = ["intel-mkl"] }
```

## How to find system MKL libraries

This crate seeks system MKL libraries, e.g. installed by [apt], [yum], or official manual installer, as following manner:

- Check `${OUT_DIR}` where previous build has downloaded
- Seek using [pkg-config] crate
  - `${PKG_CONFIG_PATH}` has to be set correctly. It may not be set by default in usual install.
  - You can confirm it by checking the following command returns error.
    ```
    pkg-config --libs mkl-dynamic-lp64-iomp
    ```
- (experimental) Seek `${XDG_DATA_HOME}/intel-mkl-tool`
- Seek a directory set by `${MKLROOT}` environment variable
- Seek default installation path
  - `/opt/intel/mkl` for Linux
  - `C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows` for Windows

If not found any MKL library and `download` feature is ON, this crate will download archive from AWS S3 `rust-intel-mkl` bucket.

[apt]: https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html
[yum]: https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-yum-repo.html
[pkg-config]: https://github.com/rust-lang/pkg-config-rs

## License
MKL is distributed under the Intel Simplified Software License for Intel(R) Math Kernel Library, See [License.txt](License.txt).
Some wrapper codes are licensed by MIT License (see the header of each file).
