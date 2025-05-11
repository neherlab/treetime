use clap::builder::styling;

pub fn styles() -> styling::Styles {
  styling::Styles::styled()
    .header(styling::AnsiColor::Green.on_default() | styling::Effects::BOLD)
    .usage(styling::AnsiColor::Green.on_default() | styling::Effects::BOLD)
    .literal(styling::AnsiColor::Blue.on_default() | styling::Effects::BOLD)
    .placeholder(styling::AnsiColor::Cyan.on_default())
}
