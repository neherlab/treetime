use app_commands::results::homoplasy::HomoplasyStatsFile;
use eyre::{Report, eyre};
use serde::de::{DeserializeOwned, IgnoredAny};
use serde_json::{Deserializer, Value};
use std::env;
use std::fs::{self, File};
use std::hint::black_box;
use std::io::{BufReader, Read};
use std::thread;
use std::time::{Duration, Instant};
use treetime_io::auspice_types::AuspiceTree;
use treetime_utils::io::compression::Decompressor;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::io::json::{json_read, json_read_file, json_read_slice};

fn main() -> Result<(), Report> {
  let args = env::args().skip(1).collect::<Vec<_>>();
  let [variants, targets, repeats, files @ ..] = args.as_slice() else {
    return Err(eyre!("usage: json_bench <variants,..> <targets,..> <repeats> <files..>"));
  };
  let repeats: usize = repeats.parse()?;
  let stack_size: usize = env::var("JSON_BENCH_STACK").map_or(Ok(2 * 1024 * 1024), |s| s.parse())?;
  for file in files {
    if env::var_os("JSON_BENCH_NOWARM").is_none() {
      black_box(fs::read(file)?);
    }
    for target in targets.split(',') {
      for variant in variants.split(',') {
        let mut best = Duration::MAX;
        for _ in 0..repeats {
          let (v, t, f) = (variant.to_owned(), target.to_owned(), file.clone());
          let elapsed = thread::Builder::new()
            .stack_size(stack_size)
            .spawn(move || measure(&v, &t, &f).map_err(|e| format!("{e:#}")))?
            .join()
            .map_err(|_| eyre!("thread panicked"))?;
          match elapsed {
            Ok(elapsed) => best = best.min(elapsed),
            Err(error) => println!("{variant:<14} {target:<8} error: {error}"),
          }
        }
        println!("{variant:<14} {target:<8} {:>8.3} s  {file}", best.as_secs_f64());
      }
    }
  }
  Ok(())
}

fn measure(variant: &str, target: &str, file: &str) -> Result<Duration, Report> {
  match target {
    "stats" => measure_typed::<HomoplasyStatsFile>(variant, file),
    "auspice" => measure_typed::<AuspiceTree>(variant, file),
    "ignored" => measure_typed::<IgnoredAny>(variant, file),
    "value" => measure_typed::<Value>(variant, file),
    _ => Err(eyre!("unknown target {target}")),
  }
}

fn measure_typed<T: DeserializeOwned>(variant: &str, file: &str) -> Result<Duration, Report> {
  let start = Instant::now();
  let value: Option<T> = match variant {
    "current" => Some(json_read_file(file)?),
    "dyn-nostack" => Some(read_nostack(open_file_or_stdin(file)?)?),
    "buf-stack" => Some(json_read(BufReader::with_capacity(256 * 1024, File::open(file)?))?),
    "bufdecomp" => {
      let file_reader = BufReader::with_capacity(256 * 1024, File::open(file)?);
      Some(json_read(BufReader::with_capacity(256 * 1024, Decompressor::from_path(file_reader, file)?))?)
    },
    "buf-nostack" => Some(read_nostack(BufReader::with_capacity(256 * 1024, File::open(file)?))?),
    "mem-ioread" => Some(json_read(fs::read(file)?.as_slice())?),
    "slice-stack" => Some(json_read_slice(&fs::read(file)?)?),
    "slice-nostack" => Some(read_slice_nostack(&fs::read(file)?)?),
    "fix" => {
      let mut bytes = Vec::new();
      open_file_or_stdin(file)?.read_to_end(&mut bytes)?;
      Some(json_read_slice(&bytes)?)
    },
    "io-only" => {
      let mut bytes = Vec::new();
      open_file_or_stdin(file)?.read_to_end(&mut bytes)?;
      black_box(bytes);
      None
    },
    "bytes-dyn" => {
      black_box(open_file_or_stdin(file)?.bytes().try_fold(0_u64, |n, b| b.map(|b| n + u64::from(b)))?);
      None
    },
    "bytes-buf" => {
      let reader = BufReader::with_capacity(256 * 1024, File::open(file)?);
      black_box(reader.bytes().try_fold(0_u64, |n, b| b.map(|b| n + u64::from(b)))?);
      None
    },
    _ => return Err(eyre!("unknown variant {variant}")),
  };
  let elapsed = start.elapsed();
  if env::var_os("JSON_BENCH_FORGET").is_some() {
    std::mem::forget(black_box(value));
  } else {
    drop(black_box(value));
  }
  Ok(elapsed)
}

fn read_nostack<T: DeserializeOwned>(reader: impl Read) -> Result<T, Report> {
  let mut de = Deserializer::from_reader(reader);
  de.disable_recursion_limit();
  let value = T::deserialize(&mut de)?;
  de.end()?;
  Ok(value)
}

fn read_slice_nostack<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, Report> {
  let mut de = Deserializer::from_slice(bytes);
  de.disable_recursion_limit();
  let value = T::deserialize(&mut de)?;
  de.end()?;
  Ok(value)
}
