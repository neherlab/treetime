use crate::cli::treetime_cli::TreetimeArgs;
use app_commands::commands::shared::argument_error::argument_error;
use clap::{Arg, Command, CommandFactory};
use clap_lex::OsStrExt;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeSet;
use std::ffi::{OsStr, OsString};
use std::iter;
use std::num::FpCategory;
use std::ops::RangeInclusive;
use treetime_utils::o;

pub(crate) const ROOT: &str = "treetime";
pub(crate) const TIMETREE: &str = "timetree";
pub(crate) const ANCESTRAL: &str = "ancestral";
pub(crate) const CLOCK: &str = "clock";
pub(crate) const MUGRATION: &str = "mugration";
pub(crate) const HOMOPLASY: &str = "homoplasy";
pub(crate) const OPTIMIZE: &str = "optimize";
pub(crate) const PRUNE: &str = "prune";
pub(crate) const PIPELINE: &str = "pipeline";

pub(crate) const EVERY_COMMAND: &[&str] = &[
  ROOT,
  TIMETREE,
  ANCESTRAL,
  CLOCK,
  MUGRATION,
  HOMOPLASY,
  OPTIMIZE,
  PRUNE,
  PIPELINE,
  "schema",
  "completions",
  "help-markdown",
  "debug",
];

const SEQUENCE_COMMANDS: &[&str] = &[TIMETREE, ANCESTRAL, HOMOPLASY, OPTIMIZE, PRUNE];
const MODEL_COMMANDS: &[&str] = &[TIMETREE, ANCESTRAL, HOMOPLASY, OPTIMIZE];
const METADATA_COMMANDS: &[&str] = &[TIMETREE, CLOCK, MUGRATION];
const SEEDED_COMMANDS: &[&str] = &[TIMETREE, ANCESTRAL, HOMOPLASY];
const OUTPUT_COMMANDS: &[&str] = &[
  TIMETREE, ANCESTRAL, CLOCK, MUGRATION, HOMOPLASY, OPTIMIZE, PRUNE, PIPELINE,
];

const MODEL_FLAGS: &[&str] = &["--gtr", "--model"];

const CLOCK_COVARIATION_INPUT: &str = "TreeTime 0.x clock reads it only for --covariation, which TreeTime 1 removed";
const COVARIATION_REMOVED: &str =
  "TreeTime 1 fits the root-to-tip regression by ordinary least squares, without the covariation model";
const COVARIATION_REROOT: &str = "this method uses the covariation model of TreeTime 0.x, which TreeTime 1 removed";

#[rustfmt::skip]
pub(crate) const V0_SPELLINGS: &[Entry] = &[
  Entry { flag: "--aln",                  arity: Arity::One,       commands: SEQUENCE_COMMANDS, rules: &[always(Then::Rename("--alignment"))] },
  Entry { flag: "--aln",                  arity: Arity::One,       commands: &[CLOCK],          rules: &[always(Then::Drop(CLOCK_COVARIATION_INPUT))] },
  Entry { flag: "--sequence-length",      arity: Arity::One,       commands: &[CLOCK],          rules: &[always(Then::Drop(CLOCK_COVARIATION_INPUT))] },
  Entry { flag: "--dates",                arity: Arity::One,       commands: METADATA_COMMANDS, rules: &[always(Then::Rename("--metadata"))] },
  Entry { flag: "--states",               arity: Arity::One,       commands: &[MUGRATION],      rules: &[always(Then::Rename("--metadata"))] },
  Entry { flag: "--name-column",          arity: Arity::One,       commands: METADATA_COMMANDS, rules: &[always(Then::Rename("--metadata-id-columns"))] },
  Entry { flag: "--outdir",               arity: Arity::One,       commands: OUTPUT_COMMANDS,   rules: &[always(Then::Rename("--output-all"))] },
  Entry { flag: "--gtr",                  arity: Arity::One,       commands: MODEL_COMMANDS,    rules: &[always(Then::Rename("--model"))] },
  Entry { flag: "--gtr-params",           arity: Arity::OneOrMore, commands: MODEL_COMMANDS,    rules: &[always(Then::Rename("--model-params"))] },
  Entry { flag: "--rng-seed",             arity: Arity::One,       commands: SEEDED_COMMANDS,   rules: &[always(Then::Rename("--seed"))] },
  Entry { flag: "--rng-seed",             arity: Arity::One,       commands: &[CLOCK],          rules: &[always(Then::Drop(CLOCK_COVARIATION_INPUT))] },
  Entry { flag: "--rng-seed",             arity: Arity::One,       commands: &[MUGRATION],      rules: &[always(Then::Drop("mugration draws no random numbers"))] },
  Entry { flag: "--aa",                   arity: Arity::None,      commands: SEQUENCE_COMMANDS, rules: &[
    when(When::Model("JTT92"), Then::Replace(&["--alphabet", "aa-no-stop"])),
    always(Then::Replace(&["--alphabet", "aa"])),
  ] },
  Entry { flag: "--keep-overhangs",       arity: Arity::None,      commands: SEQUENCE_COMMANDS, rules: &[always(Then::Replace(&["--gap-fill", "none"]))] },
  Entry { flag: "--marginal",             arity: Arity::None,      commands: &[ANCESTRAL],      rules: &[
    when(When::Given("--method-anc"), Then::Drop("--method-anc selects the method")),
    always(Then::Replace(&["--method-anc", "marginal"])),
  ] },
  Entry { flag: "--method-anc",           arity: Arity::One,       commands: &[ANCESTRAL, HOMOPLASY], rules: &[
    on(Values::OneOf(&["probabilistic", "ml"]), Then::Replace(&["--method-anc", "marginal"])),
    on(Values::OneOf(&["fitch"]),               Then::Replace(&["--method-anc", "parsimony"])),
  ] },
  Entry { flag: "--method-anc",           arity: Arity::One,       commands: &[TIMETREE],       rules: &[
    on(Values::OneOf(&["probabilistic", "ml"]), Then::Drop("timetree reconstructs ancestral sequences by marginal likelihood")),
    on(Values::OneOf(&["parsimony", "fitch"]),  Then::Stop("timetree reconstructs ancestral sequences by marginal likelihood only")),
  ] },
  Entry { flag: "--stochastic-resolve",   arity: Arity::None,      commands: &[TIMETREE],       rules: &[always(Then::Drop("TreeTime 1 resolves polytomies stochastically"))] },
  Entry { flag: "--greedy-resolve",       arity: Arity::None,      commands: &[TIMETREE],       rules: &[always(Then::Stop("TreeTime 1 resolves polytomies stochastically only"))] },
  Entry { flag: "--reroot",               arity: Arity::OneOrMore, commands: &[TIMETREE, CLOCK], rules: &[
    on(Values::OneOf(&["best", "residual", "res"]),                         Then::Replace(&["--reroot", "least-squares"])),
    on(Values::OneOf(&["min_dev"]),                                         Then::Replace(&["--reroot", "min-dev"])),
    on(Values::OneOf(&["ML"]),                                              Then::StopWith(COVARIATION_REROOT, "use --reroot least-squares")),
    on(Values::OneOf(&["min_dev_ML"]),                                      Then::StopWith(COVARIATION_REROOT, "use --reroot min-dev")),
    on(Values::OneOf(&["least-squares", "min-dev", "oldest", "midpoint"]), Then::Keep),
    always(Then::Rename("--outgroup")),
  ] },
  Entry { flag: "--coalescent",           arity: Arity::One,       commands: &[TIMETREE],       rules: &[
    on(Values::OneOf(&["const"]), Then::Replace(&["--coalescent", "opt"])),
    on(Values::Zero,              Then::Clear { setting: "coalescent", reason: "TreeTime 0.x reads 0 as no coalescent prior" }),
  ] },
  Entry { flag: "--verbose",              arity: Arity::Level,     commands: EVERY_COMMAND,     rules: &[
    on(Values::Integers(0..=0),        Then::Replace(&["--verbosity", "error"])),
    on(Values::Integers(1..=1),        Then::Replace(&["--verbosity", "warn"])),
    on(Values::Integers(2..=3),        Then::Replace(&["--verbosity", "info"])),
    on(Values::Integers(4..=u32::MAX), Then::Replace(&["--verbosity", "debug"])),
  ] },
  Entry { flag: "--time-marginal",        arity: Arity::One,       commands: &[TIMETREE],       rules: &[
    on(Values::OneOf(&["false", "never"]), Then::Drop("TreeTime 1 writes no date intervals by default")),
    on(Values::Single,                     Then::Stop("use --confidence")),
  ] },
  Entry { flag: "--branch-length-mode",   arity: Arity::One,       commands: &[TIMETREE],       rules: &[
    on(Values::OneOf(&["joint"]), Then::StopWith("TreeTime 1 estimates branch lengths by marginal likelihood", "use marginal or auto")),
  ] },
  Entry { flag: "--covariation",          arity: Arity::None,      commands: &[TIMETREE, CLOCK], rules: &[always(Then::Stop(COVARIATION_REMOVED))] },
  Entry { flag: "--tip-slack",            arity: Arity::One,       commands: &[TIMETREE, CLOCK], rules: &[always(Then::Stop(COVARIATION_REMOVED))] },
  Entry { flag: "--clock-filter-method",  arity: Arity::One,       commands: &[TIMETREE, CLOCK], rules: &[always(Then::Unavailable)] },
  Entry { flag: "--prune-outliers",       arity: Arity::None,      commands: &[TIMETREE, CLOCK], rules: &[always(Then::Rename("--prune-filtered-tips"))] },
  Entry { flag: "--plot-tree",            arity: Arity::One,       commands: &[TIMETREE],       rules: &[always(Then::Unavailable)] },
  Entry { flag: "--plot-rtt",             arity: Arity::One,       commands: &[TIMETREE, CLOCK], rules: &[always(Then::Unavailable)] },
  Entry { flag: "--tip-labels",           arity: Arity::None,      commands: &[TIMETREE],       rules: &[always(Then::Unavailable)] },
  Entry { flag: "--no-tip-labels",        arity: Arity::None,      commands: &[TIMETREE],       rules: &[always(Then::Unavailable)] },
  Entry { flag: "--n-branches-posterior", arity: Arity::None,      commands: &[TIMETREE],       rules: &[always(Then::Unavailable)] },
  Entry { flag: "--vcf-reference",        arity: Arity::One,       commands: &[TIMETREE, ANCESTRAL, CLOCK, HOMOPLASY], rules: &[always(Then::Unavailable)] },
  Entry { flag: "--custom-gtr",           arity: Arity::One,       commands: &[TIMETREE, ANCESTRAL, HOMOPLASY], rules: &[always(Then::Unavailable)] },
];

pub(crate) fn rewrite_v0_args(argv: Vec<OsString>) -> Result<V0Rewrite, Report> {
  let mut command = TreetimeArgs::command();
  command.build();

  let mut tokens = argv.into_iter();
  let mut out = Output::default();
  out.args.extend(tokens.next());
  let tokens = tokens.collect_vec();
  let separator = tokens.iter().position(|token| token == "--").unwrap_or(tokens.len());
  let (options, rest) = tokens.split_at(separator);

  let root = Scope {
    name: ROOT,
    command: &command,
  };
  let (pieces, consumed) = read_pieces(options, &root, true);
  out.rewrite(&root, pieces);
  rewrite_command(&command, &options[consumed..], &mut out);
  out.args.extend_from_slice(rest);
  out.finish()
}

#[derive(Debug, PartialEq, Eq)]
pub(crate) struct V0Rewrite {
  pub args: Vec<OsString>,
  pub messages: Vec<String>,
  pub cleared: BTreeSet<&'static str>,
}

pub(crate) struct Entry {
  pub flag: &'static str,
  pub arity: Arity,
  pub commands: &'static [&'static str],
  pub rules: &'static [Rule],
}

impl Entry {
  fn take(&self) -> Take {
    Take {
      signed_zero: self.rules.iter().any(|rule| matches!(rule.values, Values::Zero)),
      ..self.arity.take()
    }
  }
}

#[derive(Clone, Copy)]
pub(crate) enum Arity {
  None,
  One,
  OneOrMore,
  Level,
}

impl Arity {
  const fn min(self) -> usize {
    match self {
      Self::None | Self::Level => 0,
      Self::One | Self::OneOrMore => 1,
    }
  }

  const fn take(self) -> Take {
    let max = match self {
      Self::None => 0,
      Self::One | Self::Level => 1,
      Self::OneOrMore => usize::MAX,
    };
    Take {
      max,
      hyphen_values: false,
      negative_numbers: false,
      signed_zero: false,
      integers_only: matches!(self, Self::Level),
    }
  }
}

pub(crate) struct Rule {
  pub when: When,
  pub values: Values,
  pub then: Then,
}

pub(crate) enum When {
  Always,
  Given(&'static str),
  Model(&'static str),
}

pub(crate) enum Values {
  Any,
  Single,
  OneOf(&'static [&'static str]),
  Zero,
  Integers(RangeInclusive<u32>),
}

impl Values {
  fn matches(&self, values: &[OsString]) -> bool {
    match (self, values) {
      (Self::Any, _) => true,
      (_, [value]) => value.to_str().is_some_and(|value| self.matches_one(value)),
      _ => false,
    }
  }

  fn matches_one(&self, value: &str) -> bool {
    match self {
      Self::Any | Self::Single => true,
      Self::OneOf(names) => names.contains(&value),
      Self::Zero => value
        .parse::<f64>()
        .is_ok_and(|number| number.classify() == FpCategory::Zero),
      Self::Integers(range) => value.parse::<u32>().is_ok_and(|number| range.contains(&number)),
    }
  }
}

#[derive(Clone, Copy)]
pub(crate) enum Then {
  Keep,
  Rename(&'static str),
  Replace(&'static [&'static str]),
  Drop(&'static str),
  Clear {
    setting: &'static str,
    reason: &'static str,
  },
  Stop(&'static str),
  StopWith(&'static str, &'static str),
  Unavailable,
}

const fn always(then: Then) -> Rule {
  when(When::Always, then)
}

const fn when(when: When, then: Then) -> Rule {
  Rule {
    when,
    values: Values::Any,
    then,
  }
}

const fn on(values: Values, then: Then) -> Rule {
  Rule {
    when: When::Always,
    values,
    then,
  }
}

fn rewrite_command(command: &Command, tokens: &[OsString], out: &mut Output) {
  let Some((first, after)) = tokens.split_first() else {
    return;
  };
  let name = first.to_str();
  if let Some(subcommand) = name.and_then(|name| command.find_subcommand(name)) {
    out.args.push(first.clone());
    let scope = Scope {
      name: subcommand.get_name(),
      command: subcommand,
    };
    let (pieces, _) = read_pieces(after, &scope, false);
    out.rewrite(&scope, pieces);
  } else if name == Some("version") {
    out.args.push(OsString::from("--version"));
    out.args.extend_from_slice(after);
    out
      .messages
      .push(o!("TreeTime 0.x command 'version' read as '--version'"));
  } else if name == Some("arg") {
    out
      .errors
      .push(o!("the TreeTime 0.x command 'arg' is not available in TreeTime 1"));
  } else if is_flag_like(first) {
    out.args.push(OsString::from(TIMETREE));
    out.messages.push(format!(
      "TreeTime 0.x command line without a command read as '{TIMETREE}'"
    ));
    if let Some(timetree) = command.find_subcommand(TIMETREE) {
      let scope = Scope {
        name: TIMETREE,
        command: timetree,
      };
      let (pieces, _) = read_pieces(tokens, &scope, false);
      out.rewrite(&scope, pieces);
    }
  } else {
    out.args.extend_from_slice(tokens);
  }
}

fn read_pieces(tokens: &[OsString], scope: &Scope<'_>, stop_at_command: bool) -> (Vec<Piece>, usize) {
  let mut pieces = vec![];
  let mut at = 0;
  while let Some(token) = tokens.get(at) {
    let following = &tokens[at + 1..];
    let piece = match classify(token) {
      Token::Long { flag, value } => scope.long_take(flag).map(|take| {
        let values = value.map_or_else(|| take_values(following, take), |value| vec![value.to_owned()]);
        let used = if value.is_some() { 0 } else { values.len() };
        let piece = Piece::Option(Occurrence {
          flag: flag.to_owned(),
          values,
          attached: value.is_some(),
          tokens: tokens[at..=at + used].to_vec(),
        });
        (piece, 1 + used)
      }),
      Token::Short(shorts) => scope.short_take(shorts).map(|shorts| {
        let used = match shorts {
          Shorts::Complete => 0,
          Shorts::TakesValues(take) => take_values(following, take).len(),
        };
        (Piece::Plain(tokens[at..=at + used].to_vec()), 1 + used)
      }),
      Token::Value => None,
    };
    match piece {
      Some((piece, used)) => {
        pieces.push(piece);
        at += used;
      },
      None if stop_at_command => break,
      None => {
        pieces.push(Piece::Plain(vec![token.clone()]));
        at += 1;
      },
    }
  }
  (pieces, at)
}

fn classify(token: &OsStr) -> Token<'_> {
  if token.strip_prefix("--").is_some_and(|long| !long.is_empty()) {
    let (flag, value) = token
      .split_once("=")
      .map_or((token, None), |(flag, value)| (flag, Some(value)));
    return flag.to_str().map_or(Token::Value, |flag| Token::Long { flag, value });
  }
  match token.to_str().and_then(|text| text.strip_prefix('-')) {
    Some(shorts) if !shorts.is_empty() => Token::Short(shorts),
    _ => Token::Value,
  }
}

fn take_values(tokens: &[OsString], take: Take) -> Vec<OsString> {
  tokens
    .iter()
    .take(take.max)
    .take_while(|token| take.accepts(token))
    .cloned()
    .collect()
}

fn is_flag_like(token: &OsStr) -> bool {
  let bytes = token.as_encoded_bytes();
  bytes.len() > 1 && bytes.starts_with(b"-")
}

fn show(tokens: &[OsString]) -> String {
  tokens.iter().map(|token| token.to_string_lossy()).join(" ")
}

fn spelling(occurrence: &Occurrence, rule: &Rule) -> String {
  match (&rule.values, occurrence.values.as_slice()) {
    (Values::Any, _) | (_, []) => occurrence.flag.clone(),
    (_, [value]) => format!("{}={}", occurrence.flag, value.to_string_lossy()),
    (_, values) => format!("{} {}", occurrence.flag, show(values)),
  }
}

#[derive(Default)]
struct Output {
  args: Vec<OsString>,
  messages: Vec<String>,
  cleared: BTreeSet<&'static str>,
  errors: Vec<String>,
}

impl Output {
  fn rewrite(&mut self, scope: &Scope<'_>, pieces: Vec<Piece>) {
    let line = CommandLine::of(&pieces);
    let resolved = pieces
      .iter()
      .map(|piece| match piece {
        Piece::Option(occurrence) => scope.resolve(occurrence, &line),
        Piece::Plain(_) => None,
      })
      .collect_vec();
    let last_clearing = last_clearing_occurrences(&pieces, &resolved);
    for (index, (piece, rule)) in pieces.into_iter().zip(resolved).enumerate() {
      let occurrence = match piece {
        Piece::Plain(tokens) => {
          self.args.extend(tokens);
          continue;
        },
        Piece::Option(occurrence) => occurrence,
      };
      let overridden = last_clearing
        .iter()
        .any(|(flag, last)| *flag == occurrence.flag && index < *last);
      match rule {
        _ if overridden => self.messages.push(format!(
          "TreeTime 0.x argument '{}' ignored: a later {} overrides it",
          show(&occurrence.tokens),
          occurrence.flag
        )),
        Some(rule) => self.apply(occurrence, rule),
        None => self.args.extend(occurrence.tokens),
      }
    }
  }

  fn apply(&mut self, occurrence: Occurrence, rule: &Rule) {
    let original = show(&occurrence.tokens);
    match rule.then {
      Then::Keep => self.args.extend(occurrence.tokens),
      Then::Rename(target) => {
        let tokens = if occurrence.attached {
          occurrence
            .values
            .into_iter()
            .map(|value| {
              let mut token = OsString::from(format!("{target}="));
              token.push(value);
              token
            })
            .collect_vec()
        } else {
          iter::once(OsString::from(target))
            .chain(occurrence.values)
            .collect_vec()
        };
        self.messages.push(format!(
          "TreeTime 0.x argument '{original}' read as '{}'",
          show(&tokens)
        ));
        self.args.extend(tokens);
      },
      Then::Replace(targets) => {
        let tokens = targets.iter().map(OsString::from).collect_vec();
        self.messages.push(format!(
          "TreeTime 0.x argument '{original}' read as '{}'",
          show(&tokens)
        ));
        self.args.extend(tokens);
      },
      Then::Drop(reason) => {
        self
          .messages
          .push(format!("TreeTime 0.x argument '{original}' ignored: {reason}"));
      },
      Then::Clear { setting, reason } => {
        self.cleared.insert(setting);
        self.messages.push(format!(
          "TreeTime 0.x argument '{original}' clears the setting '{setting}': {reason}"
        ));
      },
      Then::Stop(advice) => {
        let spelling = spelling(&occurrence, rule);
        self
          .errors
          .push(format!("{spelling} is not available in TreeTime 1: {advice}"));
      },
      Then::StopWith(reason, advice) => {
        let spelling = spelling(&occurrence, rule);
        self
          .errors
          .push(format!("{spelling} is not available in TreeTime 1: {reason}; {advice}"));
      },
      Then::Unavailable => {
        let spelling = spelling(&occurrence, rule);
        self.errors.push(format!("{spelling} is not available in TreeTime 1"));
      },
    }
  }

  fn finish(self) -> Result<V0Rewrite, Report> {
    if self.errors.is_empty() {
      Ok(V0Rewrite {
        args: self.args,
        messages: self.messages,
        cleared: self.cleared,
      })
    } else {
      Err(argument_error(&[], self.errors.join("\n")))
    }
  }
}

fn last_clearing_occurrences(pieces: &[Piece], resolved: &[Option<&Rule>]) -> Vec<(String, usize)> {
  let clearing_flags: BTreeSet<&str> = pieces
    .iter()
    .zip(resolved)
    .filter_map(|(piece, rule)| match (piece, rule) {
      (Piece::Option(occurrence), Some(rule)) if matches!(rule.then, Then::Clear { .. }) => {
        Some(occurrence.flag.as_str())
      },
      _ => None,
    })
    .collect();
  clearing_flags
    .into_iter()
    .filter_map(|flag| {
      pieces
        .iter()
        .rposition(|piece| matches!(piece, Piece::Option(occurrence) if occurrence.flag == flag))
        .map(|last| (flag.to_owned(), last))
    })
    .collect()
}

struct Scope<'a> {
  name: &'a str,
  command: &'a Command,
}

impl Scope<'_> {
  fn entry(&self, flag: &str) -> Option<&'static Entry> {
    V0_SPELLINGS
      .iter()
      .find(|entry| entry.flag == flag && entry.commands.contains(&self.name))
  }

  fn resolve(&self, occurrence: &Occurrence, line: &CommandLine<'_>) -> Option<&'static Rule> {
    let entry = self.entry(&occurrence.flag)?;
    if occurrence.values.len() < entry.arity.min() {
      return None;
    }
    entry
      .rules
      .iter()
      .find(|rule| line.holds(&rule.when) && rule.values.matches(&occurrence.values))
  }

  fn long_take(&self, flag: &str) -> Option<Take> {
    if let Some(entry) = self.entry(flag) {
      return Some(entry.take());
    }
    let name = flag.strip_prefix("--")?;
    self
      .command
      .get_arguments()
      .find(|arg| arg.get_long() == Some(name) || arg.get_all_aliases().is_some_and(|aliases| aliases.contains(&name)))
      .map(Take::of)
  }

  fn short_take(&self, shorts: &str) -> Option<Shorts> {
    for (index, short) in shorts.char_indices() {
      let arg = self.command.get_arguments().find(|arg| {
        arg.get_short() == Some(short)
          || arg
            .get_all_short_aliases()
            .is_some_and(|aliases| aliases.contains(&short))
      })?;
      let take = Take::of(arg);
      if take.max > 0 {
        let attached = index + short.len_utf8() < shorts.len();
        return Some(if attached {
          Shorts::Complete
        } else {
          Shorts::TakesValues(take)
        });
      }
    }
    Some(Shorts::Complete)
  }
}

enum Piece {
  Plain(Vec<OsString>),
  Option(Occurrence),
}

struct Occurrence {
  flag: String,
  values: Vec<OsString>,
  attached: bool,
  tokens: Vec<OsString>,
}

enum Shorts {
  Complete,
  TakesValues(Take),
}

enum Token<'a> {
  Long { flag: &'a str, value: Option<&'a OsStr> },
  Short(&'a str),
  Value,
}

#[derive(Clone, Copy)]
struct Take {
  max: usize,
  hyphen_values: bool,
  negative_numbers: bool,
  signed_zero: bool,
  integers_only: bool,
}

impl Take {
  fn of(arg: &Arg) -> Self {
    Self {
      max: arg.get_num_args().map_or(0, |range| range.max_values()),
      hyphen_values: arg.is_allow_hyphen_values_set(),
      negative_numbers: arg.is_allow_negative_numbers_set(),
      signed_zero: false,
      integers_only: false,
    }
  }

  fn accepts(self, token: &OsStr) -> bool {
    if self.integers_only {
      return token.to_str().is_some_and(|text| text.parse::<u32>().is_ok());
    }
    !is_flag_like(token)
      || self.hyphen_values
      || (self.negative_numbers && token.to_str().is_some_and(|text| text.parse::<f64>().is_ok()))
      || (self.signed_zero && token.to_str().is_some_and(|text| Values::Zero.matches_one(text)))
  }
}

struct CommandLine<'a> {
  given: BTreeSet<&'a str>,
  models: Vec<&'a OsStr>,
}

impl<'a> CommandLine<'a> {
  fn of(pieces: &'a [Piece]) -> Self {
    let occurrences = pieces
      .iter()
      .filter_map(|piece| match piece {
        Piece::Option(occurrence) => Some(occurrence),
        Piece::Plain(_) => None,
      })
      .collect_vec();
    Self {
      given: occurrences.iter().map(|occurrence| occurrence.flag.as_str()).collect(),
      models: occurrences
        .iter()
        .filter(|occurrence| MODEL_FLAGS.contains(&occurrence.flag.as_str()))
        .flat_map(|occurrence| occurrence.values.iter().map(OsString::as_os_str))
        .collect(),
    }
  }

  fn holds(&self, when: &When) -> bool {
    match when {
      When::Always => true,
      When::Given(flag) => self.given.contains(flag),
      When::Model(name) => self
        .models
        .iter()
        .any(|model| model.to_str().is_some_and(|model| model.eq_ignore_ascii_case(name))),
    }
  }
}
