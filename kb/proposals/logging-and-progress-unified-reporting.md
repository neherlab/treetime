# Unified logging and progress reporting across packages and clients

The workspace has two independent reporting paths. Diagnostics go through the `log` crate macros, rendered by `env_logger` with the formatter in `packages/treetime-utils/src/init/global.rs`. User-facing progress goes through the `ProgressSink` trait in `packages/treetime/src/progress.rs`, with one implementation per client surface. An author writing a message has to decide which of the two to call, and a computation that reports progress has to accept a sink parameter. This proposal replaces both with one emit path built on `tracing`, where each client installs its own renderer.

## Current state

| Path                 | Emit                                   | Call sites            | Renderers                                                                                                                                                                                    |
| -------------------- | -------------------------------------- | --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Diagnostics          | `log::info!` and siblings              | 220 calls in 63 files | `env_logger`, configured from the verbosity flags in `packages/app-cli/src/cli/verbosity.rs`                                                                                                 |
| User-facing progress | `report` and `log` on a `ProgressSink` | 30 files              | `BarProgress` and `TextProgress` in `packages/app-cli/src/cli/progress.rs`, `ChannelProgress` in `packages/app-server/src/sse.rs`, `NapiProgressSink` in `packages/app-napi/src/progress.rs` |

`ProgressSink` declares `report(stage, fraction, message)`, `log(level, message)`, and `log_enabled(level)`. `NoopProgress` serves callers that report nothing.

## Motivation

- **Two paths for one concern**: the same function may need to tell the user what stage it reached and tell a developer why a value looks wrong. Those are different macros with different plumbing, and the choice is made per call site with no rule behind it
- **The sink travels through signatures**: every function on a path that reports progress carries a sink parameter, including functions that only pass it down. The parameter exists for transport, not for the computation
- **No attribution when output interleaves**: the server can run several analyses at once, and passes run in parallel through `rayon`. A diagnostic line from inside a computation cannot say which run, dataset, or partition produced it
- **Timing needs its own code**: measuring how long a pass took requires a timer written by hand at each place of interest
- **Progress reaches the interface as a sentence**: a progress event carries a stage name, a fraction, and a formatted message. Values a view could plot, such as the iteration number or the current likelihood, are either absent or buried in the text

## Examples

Each example opens with what it shows, then gives the code and the output it produces now, then the code and the output after the change, then names the packages that produce the result. Outputs are those of the current formatter in `packages/treetime-utils/src/init/global.rs`, the bar template in `packages/app-cli/src/cli/progress.rs`, and the serializers in `packages/app-server/src/sse.rs` and `packages/app-napi/src/progress.rs`.

### Diagnostic call site

A message meant for a developer, emitted from inside the root search.

Code now:

```rust
debug!("Found better node {improvements}: chi-squared improved from {best_chisq:.6e} to {tmp_chisq:.6e}");
```

Output now:

```
2026-09-19 14:03:11.482 [D] Found better node 3: chi-squared improved from 1.420000e-2 to 1.310000e-2
```

Code proposed:

```rust
debug!(improvements, from = best_chisq, to = tmp_chisq, "Found better node");
```

Output proposed:

```
2026-09-19 14:03:11.482 [D] find_best_root{run=a3f1 method=Rtt}: Found better node improvements=3 from=0.0142 to=0.0131
```

Packages: `log` with `env_logger` produce the first line. `tracing` produces the second one, rendered by the `fmt` layer of `tracing-subscriber`. The values become typed fields, the message stays constant, and the enclosing span appears before it. The layer decides how to format a field, so numeric formatting leaves the call site.

### Progress on the command line

The same stage report as today, shown as a bar in the terminal.

Code now:

```rust
pub fn reconstruct(graph: &Graph, cancel: &dyn Cancel, progress: &dyn ProgressSink) -> Result<(), Report> {
  progress.report("ancestral", 0.42, "Reconstructing internal nodes");
```

Output now:

```
⠋ [============>                 ] 42% ancestral: Reconstructing internal nodes
```

Code proposed:

```rust
pub fn reconstruct(graph: &Graph, cancel: &dyn Cancel) -> Result<(), Report> {
  info!(target: "progress", stage = "ancestral", fraction = 0.42, "Reconstructing internal nodes");
```

Output proposed:

```
⠋ [============>                 ] 42% ancestral: Reconstructing internal nodes
```

Packages: `indicatif` draws the bar in both versions. Today `BarProgress` calls it directly from the sink. In the proposed version `tracing-indicatif` owns the bar and attaches it to a span. Its `IndicatifSpanExt` trait provides `pb_set_style`, `pb_set_length`, `pb_set_position`, and `pb_set_message`, so the current model of a fraction over 1000 steps carries over without change ([span_ext.rs](https://github.com/emersonford/tracing-indicatif/blob/1853593473485272c5e7c01fa1d449beedc25f00/src/span_ext.rs)). Two details come with it: `pb_set_length` must be called before the first position update, and the `fmt` layer has to write through `indicatif_layer.get_stderr_writer()` so that diagnostics do not overwrite the bar. That writer replaces the `bar.suspend` call in `BarProgress::log`.

The displayed output is the same. The progress parameter leaves the signature; the cancellation parameter stays, because `Cancel` is a separate trait with its own purpose.

### Progress on the server stream

The same stage report, delivered to a browser or desktop client over the existing stream.

Code now:

```rust
impl ProgressSink for ChannelProgress {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    drop(self.tx.send(SinkEvent::Progress(ProgressEvent {
      stage: stage.to_owned(), fraction, message: message.to_owned(),
    })));
  }
}
```

Output now, on the wire:

```
event: progress
data: {"stage":"ancestral","fraction":0.42,"message":"Reconstructing internal nodes"}
```

Code proposed:

```rust
impl<S: Subscriber> Layer<S> for SseLayer {
  fn on_event(&self, event: &Event<'_>, ctx: Context<'_, S>) {
    let Some(progress) = ProgressEvent::from_event(event) else { return };
    let Some(run) = ctx.event_scope(event).and_then(run_id_of) else { return };
    self.route(run, SinkEvent::Progress(progress));
  }
}
```

Output proposed, on the wire:

```
event: progress
data: {"stage":"ancestral","fraction":0.42,"message":"Reconstructing internal nodes"}
```

Packages: the transport is unchanged in both versions, namely `tokio` channels, `async-stream`, and the server-sent event support of `axum`. What changes is the producer. Today `ChannelProgress` implements `ProgressSink`; in the proposed version an implementation of the `Layer` trait from `tracing-subscriber` builds the same `ProgressEvent`. The payload and the event name stay the same, so clients need no change. The routing changes: the current sink belongs to one response and needs no identifier, while the layer reads the run identifier from the span and selects the matching response channel.

### Progress on the native bridge

The same stage report, delivered to the desktop application through the native module.

Code now:

```rust
impl ProgressSink for NapiProgressSink {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    self.send_event(&NapiEvent::Progress(ProgressEvent {
      stage: stage.to_owned(), fraction, message: message.to_owned(),
    }));
  }
}
```

Output now, the string passed to the JavaScript callback:

```json
{ "type": "progress", "data": { "stage": "ancestral", "fraction": 0.42, "message": "Reconstructing internal nodes" } }
```

Code proposed:

```rust
impl<S: Subscriber> Layer<S> for NapiLayer {
  fn on_event(&self, event: &Event<'_>, _ctx: Context<'_, S>) {
    let Some(progress) = ProgressEvent::from_event(event) else { return };
    self.send_event(&NapiEvent::Progress(progress));
  }
}
```

Output proposed, the string passed to the JavaScript callback:

```json
{ "type": "progress", "data": { "stage": "ancestral", "fraction": 0.42, "message": "Reconstructing internal nodes" } }
```

Packages: the `ThreadsafeFunction` of `napi` carries the message to JavaScript and `serde_json` serializes it, in both versions. Today `NapiProgressSink` implements `ProgressSink`; in the proposed version a `Layer` from `tracing-subscriber` feeds the same call. The desktop and web code receives the same strings as before.

### Attribution under concurrency

Two analyses running at once, and the question of which one wrote a given line.

Code now, in the server:

```rust
pub fn find_best_root(graph: &Graph, params: &FindRootParams, sink: &dyn ProgressSink) -> Result<FindRootResult, Report> {
  debug!("Found better node {improvements}: chi-squared improved from {best_chisq:.6e} to {tmp_chisq:.6e}");
```

Output now, on the server's own error stream:

```
2026-09-19 14:03:11.482 [D] Found better node 3: chi-squared improved from 1.420000e-2 to 1.310000e-2
2026-09-19 14:03:11.488 [D] Found better node 1: chi-squared improved from 9.800000e-3 to 9.110000e-3
```

Neither line says which analysis produced it.

Code proposed:

```rust
#[instrument(skip(graph), fields(run = %run_id, method = ?params.method))]
pub fn find_best_root(graph: &Graph, params: &FindRootParams) -> Result<FindRootResult, Report> {
  debug!(improvements, from = best_chisq, to = tmp_chisq, "Found better node");
```

Output proposed:

```
2026-09-19 14:03:11.482 [D] find_best_root{run=a3f1 method=Rtt}: Found better node improvements=3 from=0.0142 to=0.0131
2026-09-19 14:03:11.488 [D] find_best_root{run=c8e0 method=LeastSquares}: Found better node improvements=1 from=0.0098 to=0.00911
```

Packages: the `instrument` attribute macro of `tracing` opens the span and records its fields, and the `fmt` layer of `tracing-subscriber` prints them before the message. `rayon` schedules the parallel work and does not carry the span across threads, so each worker closure has to enter it. The span also records its own duration, so per-pass timing comes from the same annotation.

### Error reports

An error raised deep in a command, and what the report tells the reader about where it happened.

Code now, in `packages/treetime-utils/src/init/global.rs`:

```rust
color_eyre::config::HookBuilder::default().install()?;
```

Output now:

```
Error: failed to read alignment

Location:
   packages/treetime-io/src/fasta/read.rs:88
```

Code proposed:

```rust
color_eyre::config::HookBuilder::default().install()?;
tracing_subscriber::registry().with(fmt_layer).with(ErrorLayer::default()).init();
```

Output proposed:

```
Error: failed to read alignment

Location:
   packages/treetime-io/src/fasta/read.rs:88

SPANTRACE:
   0: treetime::commands::ancestral with run=a3f1 dataset=flu/h3n2/20
      at packages/treetime/src/commands/ancestral.rs:41
```

Packages: `color-eyre` formats the report in both versions. The span section comes from the `ErrorLayer` of `tracing-error`, which `color-eyre` renders when the layer is installed. The exact layout is that crate's default.

### Structured progress for the interface

A report from an iterative stage, and what the desktop or web view can do with it.

Code now:

```rust
progress.report("Marginal reconstruction", 0.4, "");
```

Output now, the payload the view receives:

```json
{ "stage": "Marginal reconstruction", "fraction": 0.4, "message": "" }
```

The view can draw a bar and a label. The iteration number and the likelihood are not in the payload, and a message string would have to be parsed to recover them.

Code proposed:

```rust
info!(target: "progress", stage = "marginal", fraction = 0.4, iteration = 3, log_lh = -12043.18, "Marginal reconstruction");
```

Output proposed, the payload the view receives:

```json
{ "stage": "marginal", "fraction": 0.4, "message": "Marginal reconstruction", "fields": { "iteration": 3, "log_lh": -12043.18 } }
```

The view can plot the likelihood against the iteration while the run proceeds, and show the iteration count next to the bar.

Packages: the field values are recorded by `tracing` and read back by the layer through the visitor API of `tracing-subscriber`. This one needs a schema change: `ProgressEvent` in `packages/treetime-schema/src/progress.rs` gains a field map. The addition is compatible with the current payload, so a client that ignores the new member behaves as before. The field set is limited to names declared at the call sites, as described under Disadvantages.

### Per-stage timings

How long each stage of a command took, without writing timers.

Code now:

```rust
let started = Instant::now();
let result = reconstruct(graph, cancel, progress)?;
info!("Marginal reconstruction took {:?}", started.elapsed());
```

Output now:

```
2026-09-19 14:03:19.902 [I] Marginal reconstruction took 8.41s
```

Each stage that wants a timing needs its own pair of statements, and a stage without them reports nothing.

Code proposed:

```rust
#[instrument(skip(graph, cancel), fields(stage = "marginal"))]
pub fn reconstruct(graph: &Graph, cancel: &dyn Cancel) -> Result<(), Report> {
```

Output proposed, with the subscriber configured to print span closings:

```
2026-09-19 14:03:19.902 [I] reconstruct{run=a3f1 stage=marginal}: close time.busy=8.41s time.idle=12.4ms
```

Every instrumented stage reports its duration, and the same numbers reach a layer, so a run summary or a benchmark comparison can consume them.

Packages: the duration is measured by `tracing` on the span itself, and printed by the `fmt` layer of `tracing-subscriber` when span closings are enabled.

### Filtering one run or one stage

Selecting which output to see while a run is in progress.

Code now, the filter accepted by `env_logger`:

```
RUST_LOG=treetime::clock=debug treetime clock --tree=... --dates=...
```

Output now: every debug line from the `clock` module, from every run in the process.

Code proposed, the filter accepted by `EnvFilter`:

```
RUST_LOG=[find_best_root{method=Rtt}]=debug treetime clock --tree=... --dates=...
```

Output proposed: debug lines only from root searches that used the `Rtt` method, with lines from other methods and other runs suppressed.

Packages: the syntax for selecting on a span and its field values comes from the `EnvFilter` of `tracing-subscriber`. Module-level selection works the same way in both versions.

### Capturing events in a test

Asserting that a computation reported a warning.

Code now:

```rust
struct RecordingSink(Mutex<Vec<String>>);
impl ProgressSink for RecordingSink { /* report, log, log_enabled */ }

let sink = RecordingSink::default();
reconstruct(&graph, &NoCancel, &sink)?;
assert!(sink.0.lock().iter().any(|m| m.contains("ambiguous state")));
```

Output now: the test passes or fails, and the assertion matches on a formatted sentence, so a wording change breaks it.

Code proposed:

```rust
let events = capture_events(|| reconstruct(&graph, &NoCancel))?;
assert!(events.iter().any(|e| e.level == Level::WARN && e.field("state") == Some("ambiguous")));
```

Output proposed: the test passes or fails, and the assertion matches on a field value, so the wording of the message can change without touching the test.

Packages: the capture is a small layer on `tracing-subscriber` installed while the closure runs. It replaces the sink implementation that each test would otherwise write. A scoped capture covers one thread, so a computation that distributes work through `rayon` needs a shared collector, as described under Disadvantages.

## Wiring

The subscriber replaces `env_logger` where it is initialized today. The verbosity flags keep their shape, because `tracing` parses the same six level names and offers the same six variants, so only the imported filter type and the initialization call change:

```rust
use tracing::level_filters::LevelFilter;
tracing_subscriber::fmt().with_max_level(filter).with_writer(std::io::stderr).init();
```

The rest of the wiring:

- **Destination**: diagnostics continue to go to stderr
- **Audience marker**: the `progress` target marks events meant for the user, and events without it are diagnostics. A target was chosen over a field because a per-layer filter can select a target with a static check, while a field selector has to inspect the values of every event. The per-layer filters are what keep the progress display working under `--quiet`, so the cheaper selector is the one that supports them
- **Selection**: a renderer selects on the target, so one call site serves both audiences without the author choosing one
- **Bar ownership**: the command line bar belongs to a span, so each command opens one span that lives for the run and carries the bar
- **Dependency output**: `tracing-subscriber` captures records emitted through the `log` crate, so crates that use `log` reach the same output without change

## Disadvantages

Each item states what goes wrong, shows it, and names the way to prevent it.

### One filter governs two audiences

Verbosity today controls `env_logger` only. `ProgressSink::report` is not filtered, so the bar appears whatever the verbosity. After unification a filter sits in front of every layer.

```
treetime ancestral --quiet --tree=... --aln=...
```

Today the bar still appears and diagnostics are suppressed. With a single filter, the progress event is dropped before the bar layer sees it, and the user loses the display by asking for fewer messages.

Prevention: give each layer its own filter, so the user-facing renderers select on the `progress` target and ignore the verbosity level. This has to be part of the design, because the failure is silent.

### Routing moves from object identity to a field

The server creates a channel and a `ChannelProgress` per response, so a run's events reach the right client because the sink belongs to that response. No identifier is needed:

```rust
let (tx, rx) = mpsc::unbounded_channel::<SinkEvent>();
let progress = ChannelProgress::new(tx);
```

A single subscriber sees events from every run, so the layer has to read the run identifier from the span and select the channel:

```rust
let Some(run) = ctx.event_scope(event).and_then(run_id_of) else { return };
self.route(run, SinkEvent::Progress(progress));
```

An event emitted outside any run span has no identifier, and an event carrying the wrong one reaches another client.

Prevention: open the run span at the request boundary before any work starts, and drop events with no run identifier rather than broadcasting them.

### Parallel workers start without the span

A `rayon` worker thread has no current span:

```rust
nodes.par_iter().for_each(|node| {
  debug!(node = %node.name, "Refining branch length");
});
```

Output loses the run and stage fields for exactly the parallel work that produced the interleaving:

```
2026-09-19 14:03:12.114 [D] Refining branch length node=NODE_0000012
```

The fix is mechanical but has to be applied at every parallel boundary:

```rust
let span = Span::current();
nodes.par_iter().for_each(|node| {
  let _enter = span.enter();
  debug!(node = %node.name, "Refining branch length");
});
```

Prevention: a lint or a review checklist for `par_iter`, `par_bridge`, and `scope` call sites, since a missing entry produces output that looks correct until two runs overlap.

### Field names are fixed at the call site

Fields are declared by name where the event is written:

```rust
info!(target: "progress", stage = "marginal", iteration = 3, log_lh = -12043.18, "Marginal reconstruction");
```

A value that a view wants but that no call site declares cannot appear in the payload, and a field present only in some iterations has to be declared empty and filled in later:

```rust
let span = info_span!("marginal", iteration = field::Empty);
span.record("iteration", 3);
```

Prevention: treat the field set as part of the interface contract, listed next to the `ProgressEvent` schema, rather than as free-form data.

### Test capture covers one thread

A scoped subscriber applies to the thread that installs it:

```rust
let events = capture_events(|| reconstruct(&graph, &NoCancel))?;
```

If `reconstruct` distributes work through `rayon`, the events emitted on worker threads are not captured, and the assertion passes or fails for the wrong reason.

Prevention: capture through a collector shared across threads, or restrict the capture helper to computations that stay on one thread and say so where it is defined.

### An attribute without `skip` formats its arguments

```rust
#[instrument]
pub fn reconstruct(graph: &Graph, cancel: &dyn Cancel) -> Result<(), Report> {
```

Every call formats the whole graph through `Debug` to record it as a span field. The cost appears at runtime and nothing reports it as a mistake.

Prevention: `skip` every argument that is not a small scalar or a name, and state that rule where instrumentation is described.

### Hot paths must stay free of spans

A span per tree node or per alignment site costs measurable time. Instrumentation belongs at command, pass, and request level, and a timing comparison against the current code is the check that this held.

### Smaller items

- **A failing layer runs on the computation thread**: a panic inside `on_event` takes down the worker, while a failed channel send is dropped today
- **Silence has three causes instead of one**: a missing subscriber, a filter, or layer ordering. Diagnosing takes longer than with `env_logger`

## Alternatives

- **Convert diagnostics only**: replace `log` with `tracing` and keep `ProgressSink` for user-facing progress. Smaller change, and routing stays with the response object, but the two paths and the sink parameter remain
- **Extend `ProgressSink` with structured fields and a run identifier**: keeps the current transport and adds attribution to it. This reimplements span propagation, field recording, and filtering inside the project
- **Keep both paths unchanged**: no work, and the attribution problem stays. Acceptable only if the server never runs analyses concurrently

## Migration path

1. **Add the dependencies**: `tracing`, `tracing-subscriber`, and `tracing-error`, and initialize the subscriber where `env_logger` is initialized today
2. **Convert the diagnostic call sites**: all 220 of them, turning values formatted into the message into fields
3. **Add the spans**: at command entry points and pass boundaries, carrying the run identifier, the dataset, and the stage
4. **Cover the parallel boundaries**: enter the captured span inside every `rayon` worker closure
5. **Convert the progress call sites**: `ProgressSink::report` and `ProgressSink::log` become events on the `progress` target
6. **Replace the sinks with layers**: one layer per client, then remove the trait, the macros, and the sink parameters
7. **Remove the old dependencies**: `log` and `env_logger` leave the workspace

Steps 3 and 4 belong together. A span added without covering the parallel boundaries below it reports attribution that is missing from exactly the output that needs it.

## Verification

- **Routing**: two concurrent server runs produce two event streams, each carrying only its own run
- **Client output**: the command line bar, the stream payloads, and the bridge callback strings match the outputs shown above, byte for byte where they are serialized
- **Verbosity**: the flags select the same levels as before for each command
- **Cost**: a representative pass shows no measurable slowdown, confirming that no span sits in a hot loop

## Open questions

- **Subscriber ownership on the bridge**: whether the native module may install a process-wide subscriber, given that the host process belongs to the desktop application. Nothing else in that process installs one today, so the question is one of ownership rather than conflict. If the answer is yes, the module installs it once from its initializer so that repeated calls do nothing

## Evidence still to gather

A spike wiring one command end to end would settle the routing and filter design before the migration starts. It covers one command, one span at the entry point, the bar layer, the stream layer, and two concurrent runs against the server. It answers whether routing by span field delivers each run to the right client, whether entering the span in `rayon` workers restores the fields, whether per-layer filters keep the bar under `--quiet`, and what a span closing costs on a representative pass. The spike is throwaway code and commits nothing.
