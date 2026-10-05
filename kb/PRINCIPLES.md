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
    pub seed: u64,
  }

  pub fn run(input: &Input, params: &Params, sink: &mut dyn Sink) -> Result<Output, Error> {
    let c = plan(input, params);
    emit_items(&input.right, &c, params.seed, sink)?;
    let e = E::draw(&c, params.seed);
    Ok(Output::assemble(&input.right, &c, &e))
  }

  pub trait Sink {
    fn emit(&mut self, item: Item) -> Result<(), Error>;
  }

  fn plan(input: &Input, params: &Params) -> C {
    let a = A::build(&input.left, params.tolerance);
    let b = step_b(&input.right, &a);
    let d = collect_d(&input.parts, &b, params.tolerance);
    step_c(&a.secondary, &b, &d, params)
  }

  fn emit_items(right: &Right, c: &C, seed: u64, sink: &mut dyn Sink) -> Result<(), Error> {
    for item in Item::stream(right, c, seed) {
      sink.emit(item)?;
    }
    Ok(())
  }

  fn step_b(right: &Right, a: &A) -> B {
    match a.mode {
      Mode::Fast => B::fast(right),
      Mode::Full => B::build(right, &a.primary),
    }
  }

  fn collect_d(parts: &[Part], b: &B, tolerance: f64) -> Vec<D> {
    parts
      .iter()
      .filter(|part| part.selected(tolerance))
      .map(|part| D::build(part, b))
      .collect()
  }

  fn step_c(secondary: &Secondary, b: &B, d: &[D], params: &Params) -> C {
    let seed = Seed::build(secondary, d);
    let weights = Weights::build(b, d);

    let mut state = State::init(&seed, &weights);
    for round in 0..params.iterations {
      let proposal = state.advance(&weights, round);
      let scored = Scored::eval(&proposal, b);
      state = if scored.improves(&state, params.tolerance) {
        state.accept(scored)
      } else {
        state.relax(proposal)
      };
      if state.converged(params.tolerance) {
        break;
      }
    }
    C::finalize(state, &seed)
  }

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

mod format {
  pub fn parse(raw: &[u8]) -> Result<core::Input, Error> { core::Input::parse(raw) }
  pub fn encode(output: &core::Output) -> Vec<u8> { output.encode() }
  pub fn encode_item(item: &core::Item) -> Vec<u8> { item.encode() }
}

mod validation {
  pub fn validate(params: &core::Params) -> Result<(), Error> { params.check() }
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
    pub summary: Option<PathBuf>,
    #[arg(long, default_value_t = 10)]
    pub iterations: usize,
    #[arg(long, default_value_t = 1e-6)]
    pub tolerance: f64,
    #[arg(long)]
    pub seed: Option<u64>,
  }

  impl From<&CliArgs> for core::Params {
    fn from(a: &CliArgs) -> Self {
      core::Params { iterations: a.iterations, tolerance: a.tolerance, seed: a.seed.unwrap_or_else(draw_seed) }
    }
  }

  pub fn run(args: &CliArgs) -> Result<(), Error> {
    let input = format::parse(&fs::read(&args.input)?)?;
    let params = core::Params::from(args);
    validation::validate(&params)?;
    let mut sink = FileSink { out: File::create(&args.output)? };
    let output = core::run(&input, &params, &mut sink)?;
    if let Some(path) = &args.summary {
      fs::write(path, format::encode(&output))?;
    }
    Ok(())
  }

  struct FileSink {
    out: File,
  }

  impl core::Sink for FileSink {
    fn emit(&mut self, item: core::Item) -> Result<(), Error> {
      self.out.write_all(&format::encode_item(&item))?;
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
    let input = format::parse(&req.input)?;
    validation::validate(&req.params)?;
    let mut sink = ResponseSink { body };
    let output = core::run(&input, &req.params, &mut sink)?;
    sink.body.write_chunk(&format::encode(&output))
  }

  struct ResponseSink {
    body: ResponseStream,
  }

  impl core::Sink for ResponseSink {
    fn emit(&mut self, item: core::Item) -> Result<(), Error> {
      self.body.write_chunk(&format::encode_item(&item))
    }
  }
}
```

## Principle 1: Inference is a value pipeline with no shared mutable state.

- **P1.1. Value pipeline**:
  - **P1.1.1. Inference as a chain**: inference is a chain of `step(input1, input2, ...) -> output` functions, organized hierarchically
  - **P1.1.2. Per-item data as value maps**: per-item data flows as value collections keyed by a stable identifier, not by mutating shared state
  - **P1.1.3. DAG dataflow**: data flows along a DAG: sources -> transformations -> sinks, and only that direction. Sources are read, never written back. Each stage is a pure function that consumes values and returns new values; no stage mutates its inputs, a shared object, or a neighbor's state. The dataflow is acyclic: nothing feeds back upstream, and only the sink is written. A single long-lived mutable object spanning stages is the pattern this forbids
- **P1.2. No god-objects**: no long-lived mutable object spans stages, and no object is partially initialized then completed by later mutation. State flows as values passed in and returned. The values need not mirror any earlier in-place structure; they carry only the necessary information. A step takes several inputs or one parameter record (P2.5), and may use its inputs and prior outputs fully or partially. The final sink is the output write.
- **P1.3. Structure without payload**: a shared structural topology object (graph, tree) carries only its shape, holding no per-item data payload and no data-type generics. Per-item data travels alongside it.
- **P1.4. External iteration**: a traversal returns the visit order (`Vec<NodeKey>` or an iterator), and the caller loops with `for`. Reason: `?`, `break`, borrowing and composition work in a loop, not inside a callback. Exception: a parallel scheduler that owns its work queue.
- **P1.5. Complete maps**: every per-item map covers every item it describes. A missing entry is a bug: index with `map[&key]`, never `get(..).unwrap_or_default()`.

## Principle 2: Top-level functions read as an ordered sequence of named, single-responsibility steps.

- **P2.1. Orchestrator, not a blob**: a top-level function is a short orchestrator whose body reads top-to-bottom like a table of contents. It names the steps and fixes their order; it does not inline their guts. A reader traces the whole operation from the orchestrator alone, then descends only into the one step that matters. Branches and loops are a step's internals, never the orchestrator body: an orchestrator that must choose or iterate names a selecting or iterating step, so the top level stays a flat sequence of `let x = step(...)` calls. Exception: a flat list of optional outputs (`if let Some(path) = &args.x { write_x(path, &value)? }`) or one `match` over requested outputs stays in the orchestrator, because extracting it creates a single-caller wrapper (P5.1). A step may itself be a sub-orchestrator with its own table of contents (see `step_c` in the example).
- **P2.2. One responsibility per step**: each step does one thing -- read, one transformation, or write -- never a mix. The phases stay separated and in order: inputs read first, computation in the middle, outputs written last. No write concern leaks into the compute phase, and no I/O hides inside a computational step.
- **P2.3. Correct, minimal boundaries**: each step takes exactly the inputs it needs and returns exactly what its consumers read, no more. No pass-through parameters (threaded in and handed back unused), no unused return fields (produced and never read), no god-struct passed so helpers can reach unrelated state.
- **P2.4. Pass data, not behavior**: a step receives the data it needs as maps, slices or plain structs. Use a closure only as a call-site argument to a library combinator (`map`, `sort_by_key`) or a scoped resource (`write_file_with(path, |f| ..)`). Use a trait object only for an open set with two or more production implementations chosen at run time (the sink in the example). Reason: closures and trait objects hide what a step reads, block `?` and `break`, and let per-caller differences grow unseen.
- **P2.5. Parameter records**: an immutable struct that groups the inputs of one step, built once and passed by reference, is a parameter record, not a god-object (P1.2) or a command object (P3.8). Prefer it over more than five positional parameters, or when several callers pass the same group. Group optional data that belongs together in one `Option<Group>` or an enum, so a value cannot exist without the values it depends on (mutations without their reference sequence).
- **P2.6. Linear passes**: per-item work is O(n log n) or better: no walk to the root per node, no linear search per node, no output per pair of items.

## Principle 3: Separation of Concerns

One core, served unchanged from a CLI and a web backend. The CLI and web adapters own their own I/O and never depend on each other.

- **P3.1. Pure core**: core algorithms accept input parameters, perform computations, and return results without relying on or modifying external state, I/O, CLI, or Web context
- **P3.2. CLI adapter**: the CLI layer handles CLI args and I/O (files, stdin, stdout etc.) and its interaction with the core algorithms, without embedding algo logic
- **P3.3. Web adapter**: the web backend layer handles HTTP requests, responses, and its interaction with the core algorithms, without embedding algo logic
- **P3.4. Shared is a tier, not a module**: code used by more than one adapter (parse/encode of core types, shared validation) lives in the shared tier -- below the adapters, beside the core. Split it into cohesive units named for their concern (`newick`, `auspice`, `output-paths`), one reason to change each. NEVER one catch-all crate or module; NEVER a layer or grab-bag name (`shared`, `common`, `api`, `util`, `misc`). Adapter-specific (de)serialization -- CLI argument parsing, HTTP request and response shapes -- stays in its own adapter
- **P3.5. One-way dependencies**: dependencies point one way and form no cycle. The CLI and web adapters depend on the core and on shared code; the core depends on neither adapter and holds no CLI or Web knowledge; the two adapters never depend on each other
- **P3.6. Return values; stream only large per-item results**: functions return their results by value. Only a core `run` that produces large per-item output while it computes (reconstructed sequences, optimizer trace) emits it in order through a caller-supplied sink, so the core owns no I/O and holds no full copy. The CLI writes the stream to a file, the web backend to a response stream. Output projection, graph traversal and every other step return values (P1.4, P2.4). Never pass raw `Write` or bytes into the core interface, because that leaks format and I/O into the core.
- **P3.7. No shared blob**: one shared unit per concern, format, or domain, so removing a capability touches one unit and a consumer depends only on what it uses. The no-`utils` rule applies inside the shared tier exactly as everywhere else.
- **P3.8. Adapters own the workflow, shared code owns steps**: each adapter owns its complete external interaction -- convert its arguments or request, read the inputs, build core values, call the core, then project and deliver the result. A shared unit owns one step (parse, validate, plan, project, encode, write). It never owns a whole read-run-write command or an adapter-neutral command object, because that recreates the generic command layer under a new name. A plain data struct that several writers read, such as a tree with its per-node facts, is a parameter record (P2.5), not a command object.
- **P3.9. The encoding tier encodes only**: a shared output unit turns core values into bytes and nothing more. It reads no input, parses nothing, and holds no path policy. The adapter selects the destination; the writer may open the path the adapter chose.

## Principle 4: Production code is the only consumer that shapes structure.

Production code is the sole consumer that shapes architecture and refactoring

- **P4.1. Production call sites decide**: module boundaries, interfaces, signatures, visibility, and placement derive from production callers only
- **P4.2. Tests do not count**: a symbol referenced solely by tests has zero consumers; exclude test files from fan-in, single-consumer chain, and unused-code counts
- **P4.3. Tests never justify a shape**: a test, fixture, mock, or helper is never a reason to add, keep, or shape a production unit, parameter, boundary
- **P4.4. Tests come last**: the architecture is settled from production first, and after that the tests are updated, moved, or deleted to fit the new structure

## Principle 5: Design from intent, with no accidental abstractions or legacy remnants.

- **P5.1. No accidental abstractions**: distrust incidental structure. Remove a wrapper with a single caller, a struct destructured one line after it is built, an `Option` that is always `Some` on success, and a type that duplicates an existing canonical type. Keep duplication until a stable shared concept appears; merge entangled or single-consumer code, not merely similar code. When three or more copies differ only in data, put the data in a struct and share one function, because such copies drift apart; copies that differ in logic stay separate.
- **P5.2. Design from intent**: name and shape each unit from what the operation must do, not from an incidental code layout. Prefer a canonical type or utility over a new local one. The result reads as written from scratch: no versioned names, no edit-history comments, no compatibility shims, legacy wrappers, or re-exports to a prior design.
- **P5.3. Typed values until the encoder**: never format a value to a string and parse it back; carry the typed value to the writer. Reason: a round trip through text changes values (`2020.50` becomes `2020.5`, the trait value `Nan` becomes `NaN`).

## Principle 6: Each operation exposes one uniform core contract.

Each operation in the core has one entry point and one result, named and shaped the same way across operations, so a reader locates any operation's contract without opening a file.

- **P6.1. One canonical result**: each operation owns exactly one aggregate `Output` value in the core. Every response body, file format, and wire form is a projection derived from it in an adapter or the shared tier, never a second result type kept in the core.
- **P6.2. Uniform naming**: each operation exposes the same-named `Params`, `Input`, `Output`, and `run`. `Input` holds fully parsed domain values, not raw arguments, paths, or parser records. Rename or relocate a single-consumer type to fit this shape rather than wrap it in a new one.
- **P6.3. Uniform shape, not uniform parameters** (extends P2.3): uniformity governs names and order, never the presence of a parameter. Add a sink, progress, or cancellation argument only to an operation that has that concern; a matching signature is never a reason to carry an argument the operation ignores.
- **P6.4. Signals split by direction**: cancellation is a read-only input the computation polls; progress and trace are output sinks the computation writes. They are separate parameters, never methods on one trait and never bundled into a context object. Sinks carry domain values, not encoded bytes or paths (see P3.6).
- **P6.5. Typed failures keep the cause chain**: an operation returns typed failure classes at its boundary. Classification adds meaning; it never rewrites, truncates, or flattens the underlying message chain.
