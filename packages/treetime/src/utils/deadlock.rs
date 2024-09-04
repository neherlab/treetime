use ndarray_linalg::Transpose::No;

pub fn enable_deadlock_detection() {
  #[cfg(any(debug_assertions, test, feature = "deadlock_detection"))]
  {
    use crate::utils::global_init::{HIDDEN_CRATE_NAME_PREFIXES, HIDDEN_CRATE_PATH_PREFIXES};
    use color_backtrace::{
      termcolor::{Color, ColorChoice, ColorSpec, StandardStream},
      BacktracePrinter, ColorScheme,
    };
    use color_eyre::owo_colors::OwoColorize;
    use parking_lot::deadlock;
    use std::thread;
    use std::time::Duration;

    fn cs(fg: Option<Color>, intense: bool, bold: bool) -> ColorSpec {
      let mut cs = ColorSpec::new();
      cs.set_fg(fg);
      cs.set_intense(intense);
      cs.set_bold(bold);
      cs
    }

    let backtrace_printer = BacktracePrinter::new()
      .color_scheme(ColorScheme {
        frames_omitted_msg: cs(Some(Color::Black), false, true),
        src_loc: cs(Some(Color::Green), false, false),
        src_loc_separator: cs(Some(Color::Green), false, false),
        dependency_code: cs(Some(Color::Green), false, false),
        crate_code: cs(Some(Color::Red), true, false),
        selected_src_ln: cs(Some(Color::Cyan), false, true),
        ..ColorScheme::default()
      })
      .strip_function_hash(true)
      .add_frame_filter(Box::new(|frames| {
        frames.retain(|frame| {
          let should_show_name = frame.name.as_ref().map_or(false, |name| {
            !HIDDEN_CRATE_NAME_PREFIXES
              .iter()
              .any(|&prefix| name.starts_with(prefix) || name.starts_with(&format!("<{prefix}")))
          });

          let should_show_file = !frame.filename.as_ref().map_or(false, |filename| {
            HIDDEN_CRATE_PATH_PREFIXES
              .iter()
              .any(|&prefix| filename.starts_with(prefix))
          });

          should_show_file && should_show_name
        });
      }));

    let mut color_stderr = StandardStream::stderr(ColorChoice::Always);

    #[allow(clippy::infinite_loop)]
    thread::spawn(move || loop {
      thread::sleep(Duration::from_secs(5));

      let deadlocks = deadlock::check_deadlock();
      if deadlocks.is_empty() {
        continue;
      }

      eprintln!(
        "{}",
        format!("{} deadlock(s) detected", deadlocks.len()).bright_red().bold()
      );
      for (i, threads) in deadlocks.iter().enumerate() {
        eprintln!("{}", format!("Deadlock #{i}").bright_yellow().bold());
        for t in threads {
          eprintln!("{}", format!("Thread Id {:#?}", t.thread_id()).bright_cyan());
          backtrace_printer
            .print_trace(t.backtrace(), &mut color_stderr)
            .expect("Failed to print backtrace");
        }
      }
    });
  }
}
