use std::process::ExitCode;

const INSTALL_HELP: &str = r#"TreeTime links OpenBLAS through pkg-config, which did not find it. Install OpenBLAS and pkg-config:

  Debian, Ubuntu:   apt install libopenblas-dev pkg-config
  macOS:            brew install openblas pkgconf
                    export PKG_CONFIG_PATH="$(brew --prefix openblas)/lib/pkgconfig"
  FreeBSD:          pkg install openblas pkgconf
  Windows (GNU):    pacman -S mingw-w64-ucrt-x86_64-openblas mingw-w64-ucrt-x86_64-pkgconf  (MSYS2)
  Windows (MSVC):   conda install -c conda-forge openblas pkg-config

Homebrew keeps OpenBLAS out of its main prefix, so pkg-config finds it only through PKG_CONFIG_PATH.
See the section "Host" in docs/dev/developer_guide.md for builds outside the build container."#;

fn main() -> ExitCode {
  match pkg_config::Config::new().probe("openblas") {
    Ok(_) => ExitCode::SUCCESS,
    Err(error) => {
      eprint!("{error}\n\n{INSTALL_HELP}\n");
      ExitCode::FAILURE
    },
  }
}
