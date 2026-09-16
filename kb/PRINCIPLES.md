# Principles

These are the guiding principles, kept as a reference to ensure consistency and prevent drift while implementing and cleaning up. Read each item as a standard to hold work against, not a status checklist: principles already satisfied in the code stay here as guardrails against regression, and the rest direct the work still in progress.

## Idealistic example

Abstract example which applies all principles. Real implementation might have more nuance - to be brought forward prominently, discussed and handled explicitly.

```rust
mod core {
  #[derive(Clone, Debug, Serialize, Deserialize)]
  pub struct Params {
    pub iterations: usize,
    pub tolerance: f64,
    pub seed: Option<u64>,
  }

  pub fn run(input: &Input, params: &Params, sink: &mut dyn Sink) -> Result<Output, Error> {
    let c = plan(input, params);
    emit_items(&input.right, &c, params.seed, sink)?;
    let e = step_e(&c, params.seed);
    Ok(step_f(&input.right, &c, &e))
  }

  pub trait Sink {
    fn emit(&mut self, item: Item) -> Result<(), Error>;
  }

  fn plan(input: &Input, params: &Params) -> C {
    let a = step_a(&input.left, params.tolerance);
    let b = step_b(&input.right, &a);
    let d = collect_d(&input.parts, &b, params.tolerance);
    step_c(&a.secondary, &b, &d, params)
  }

  fn emit_items(right: &Right, c: &C, seed: Option<u64>, sink: &mut dyn Sink) -> Result<(), Error> {
    for item in Item::stream(right, c, seed) {
      sink.emit(item)?;
    }
    Ok(())
  }

  fn step_a(left: &Left, tolerance: f64) -> A { A::build(left, tolerance) }

  fn step_b(right: &Right, a: &A) -> B {
    match a.mode {
      Mode::Fast => B::fast(right),
      Mode::Full => B::build(right, &a.primary),
    }
  }

  fn collect_d(parts: &[Part], b: &B, tolerance: f64) -> Vec<D> {
    let mut items = Vec::new();
    for part in parts {
      if part.selected(tolerance) {
        items.push(build_d(part, b));
      }
    }
    items
  }

  fn build_d(part: &Part, b: &B) -> D { D::build(part, b) }

  fn step_e(c: &C, seed: Option<u64>) -> E { E::draw(c, seed) }
  fn step_f(right: &Right, c: &C, e: &E) -> Output { Output::assemble(right, c, e) }

  fn step_c(secondary: &Secondary, b: &B, d: &[D], params: &Params) -> C {
    let seed = seed_c(secondary, d);
    let weights = weigh(b, d);

    let mut state = init_c(&seed, &weights);
    for round in 0..params.iterations {
      let proposal = advance_c(&state, &weights, round);
      let scored = score_c(&proposal, b);
      state = if scored.improves(&state, params.tolerance) {
        accept_c(state, scored)
      } else {
        relax_c(state, proposal)
      };
      if state.converged(params.tolerance) {
        break;
      }
    }
    finalize_c(state, &seed)
  }

  fn seed_c(secondary: &Secondary, d: &[D]) -> Seed { Seed::build(secondary, d) }
  fn weigh(b: &B, d: &[D]) -> Weights { Weights::build(b, d) }
  fn init_c(seed: &Seed, weights: &Weights) -> State { State::init(seed, weights) }
  fn advance_c(state: &State, weights: &Weights, round: usize) -> Proposal { state.advance(weights, round) }
  fn score_c(proposal: &Proposal, b: &B) -> Scored { Scored::eval(proposal, b) }
  fn accept_c(state: State, scored: Scored) -> State { state.accept(scored) }
  fn relax_c(state: State, proposal: Proposal) -> State { state.relax(proposal) }
  fn finalize_c(state: State, seed: &Seed) -> C { C::finalize(state, seed) }

  enum Mode {
    Fast,
    Full,
  }

  struct A {
    mode: Mode,
    primary: Primary,
    secondary: Secondary,
  }
}

mod shared {
  pub fn parse(raw: &[u8]) -> Result<core::Input, Error> { core::Input::parse(raw) }
  pub fn validate(params: &core::Params) -> Result<(), Error> { params.check() }
  pub fn encode(output: &core::Output) -> Vec<u8> { output.encode() }
  pub fn encode_item(item: &core::Item) -> Vec<u8> { item.encode() }
}

mod cli {
  #[derive(clap::Parser)]
  #[command(rename_all = "kebab-case")]
  pub struct CliArgs {
    #[arg(value_hint = ValueHint::FilePath)]
    pub input: PathBuf,
    #[arg(long, value_hint = ValueHint::FilePath)]
    pub output: PathBuf,
    #[arg(long, value_hint = ValueHint::FilePath)]
    pub summary: PathBuf,
    #[arg(long, default_value_t = 10)]
    pub iterations: usize,
    #[arg(long, default_value_t = 1e-6)]
    pub tolerance: f64,
    #[arg(long)]
    pub seed: Option<u64>,
  }

  impl From<&CliArgs> for core::Params {
    fn from(a: &CliArgs) -> Self {
      core::Params { iterations: a.iterations, tolerance: a.tolerance, seed: a.seed }
    }
  }

  pub fn run(args: &CliArgs) -> Result<(), Error> {
    let input = shared::parse(&fs::read(&args.input)?)?;
    let params = core::Params::from(args);
    shared::validate(&params)?;
    let mut sink = FileSink { out: File::create(&args.output)? };
    let output = core::run(&input, &params, &mut sink)?;
    fs::write(&args.summary, shared::encode(&output))?;
    Ok(())
  }

  struct FileSink {
    out: File,
  }

  impl core::Sink for FileSink {
    fn emit(&mut self, item: core::Item) -> Result<(), Error> {
      self.out.write_all(&shared::encode_item(&item))?;
      Ok(())
    }
  }
}

mod web {
  #[derive(Deserialize)]
  pub struct Request {
    pub input: Vec<u8>,
    pub params: core::Params,
  }

  pub fn handle(req: &Request, body: ResponseStream) -> Result<(), Error> {
    let input = shared::parse(&req.input)?;
    shared::validate(&req.params)?;
    let mut sink = ResponseSink { body };
    let output = core::run(&input, &req.params, &mut sink)?;
    sink.body.write_chunk(&shared::encode(&output))
  }

  struct ResponseSink {
    body: ResponseStream,
  }

  impl core::Sink for ResponseSink {
    fn emit(&mut self, item: core::Item) -> Result<(), Error> {
      self.body.write_chunk(&shared::encode_item(&item))
    }
  }
}
```

## Principle 1: Inference is a value pipeline with no shared mutable state.

- **P1.1. Value pipeline**:
  - **P1.1.1. Inference as a chain**: inference is a chain of `step(input1, input2, ...) -> output` functions, organized hierarchically
  - **P1.1.2. Per-item data as value maps**: per-item data flows as value collections keyed by a stable identifier, not by mutating shared state
  - **P1.1.3. DAG dataflow**: data flows along a DAG: sources -> transformations -> sinks, and only that direction. Sources are read, never written back. Each stage is a pure function that consumes values and returns new values; no stage mutates its inputs, a shared object, or a neighbor's state. The dataflow is acyclic: nothing feeds back upstream, and only the sink is written. A single long-lived mutable object spanning stages is the pattern this forbids
- **P1.2. No god-objects**: no long-lived mutable object spans stages, and no object is partially initialized then completed by later mutation. State flows as values passed in and returned. The values need not mirror any earlier in-place structure; they carry only the necessary information. Multiple inputs may be passed rather than one object, and a step may use its inputs and prior outputs fully or partially. The final sink is the output write.
- **P1.3. Structure without payload**: a shared structural object carries only its shape, holding no per-item data payload and no data-type generics. Per-item data travels as value maps alongside it.

## Principle 2: Top-level functions read as an ordered sequence of named, single-responsibility steps.

- **P2.1. Orchestrator, not a blob**: a top-level function is a short orchestrator whose body reads top-to-bottom like a table of contents. It names the steps and fixes their order; it does not inline their guts. A reader traces the whole operation from the orchestrator alone, then descends only into the one step that matters. Branches and loops are a step's internals, never the orchestrator body: an orchestrator that must choose or iterate names a selecting or iterating step, so the top level stays a flat sequence of `let x = step(...)` calls. A step may itself be a sub-orchestrator with its own table of contents (see `step_c` in the example).
- **P2.2. One responsibility per step**: each step does one thing -- read, one transformation, or write -- never a mix. The phases stay separated and in order: inputs read first, computation in the middle, outputs written last. No write concern leaks into the compute phase, and no I/O hides inside a computational step.
- **P2.3. Correct, minimal boundaries**: each step takes exactly the inputs it needs and returns exactly what its consumers read, no more. No pass-through parameters (threaded in and handed back unused), no unused return fields (produced and never read), no god-struct passed so helpers can reach unrelated state.

## Principle 3: Separation of Concerns

One core, served unchanged from a CLI and a web backend. The CLI and web adapters own their own I/O and never depend on each other.

- **P3.1. Pure core**: core algorithms accept input parameters, perform computations, and return results without relying on or modifying external state, I/O, CLI, or Web context
- **P3.2. CLI adapter**: the CLI layer handles CLI args and I/O (files, stdin, stdout etc.) and its interaction with the core algorithms, without embedding algo logic
- **P3.3. Web adapter**: the web backend layer handles HTTP requests, responses, and its interaction with the core algorithms, without embedding algo logic
- **P3.4. Shared code**: functionality common to more than one adapter (parsing and encoding of core types, shared validation) lives in a shared layer, written once. Adapter-specific (de)serialization, such as CLI argument parsing and HTTP request and response shapes, stays in its own adapter
- **P3.5. One-way dependencies**: dependencies point one way and form no cycle. The CLI and web adapters depend on the core and on shared code; the core depends on neither adapter and holds no CLI or Web knowledge; the two adapters never depend on each other
- **P3.6. Return small results, stream large ones**: one `run` returns the small aggregate output by value and emits the large per-item output in order through a caller-supplied sink; the core owns no I/O either way. The CLI writes the stream to a file, the web backend to a response stream. The sink can take several forms:
  - **P3.6.1. Push callback / trait sink** (recommended): matches existing code, core owns no I/O
  - **P3.6.2. Returned lazy iterator (pull)**: cleanest separation of concerns, but Rust lifetime friction
  - **P3.6.3. Channel**: good for parallel ordered emit, adds concurrency machinery
  - **P3.6.4. Async stream**: best for web streaming, risks async in core
  - **P3.6.5. Raw `Write` of bytes**: reject at the core seam (leaks format/I/O)

## Principle 4: Production code is the only consumer that shapes structure.

Production code is the sole consumer that shapes architecture and refactoring

- **P4.1. Production call sites decide**: module boundaries, seams, signatures, visibility, and placement derive from production callers only
- **P4.2. Tests do not count**: a symbol referenced solely by tests has zero consumers; exclude test files from fan-in, single-consumer chain, and unused-code counts
- **P4.3. Tests never justify a shape**: a test, fixture, mock, or helper is never a reason to add, keep, or shape a production unit, parameter, boundary
- **P4.4. Tests come last**: the architecture is settled from production first, and after that the tests are updated, moved, or deleted to fit the new structure

## Principle 5: Design from intent, with no accidental abstractions or legacy remnants.

- **P5.1. No accidental abstractions**: distrust incidental structure. Remove a wrapper with a single caller, a struct destructured one line after it is built, an `Option` that is always `Some` on success, and a type that duplicates an existing canonical type. Keep duplication until a stable shared concept appears; merge entangled or single-consumer code, not merely similar code.
- **P5.2. Design from intent**: name and shape each unit from what the operation must do, not from an incidental code layout. Prefer a canonical type or utility over a new local one. The result reads as written from scratch: no versioned names, no edit-history comments, no compatibility shims, legacy wrappers, or re-exports to a prior design.

## Principle 6: Behavior is a specified contract, changed only deliberately.

- **P6.1. Behavior is a specified contract**: the observable output is defined, not incidental; it stays stable except where a change is deliberate and stated.
- **P6.2. Observable changes are deliberate**: any change to observable behavior is surfaced and approved before it lands, never made silently.
