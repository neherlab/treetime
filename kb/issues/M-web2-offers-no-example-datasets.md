# Web2 offers no example datasets

The back end of web2 serves an empty examples folder (`/treetime/examples`, created by `Backend::new()` in `packages/app-wasm/src/backend.rs`), so the dataset catalog of the interface is empty and a user starts every run from uploaded files.

## Evidence

- The examples download of the local apps (`packages/app-commands/src/examples_fetch.rs`) uses the blocking client of `reqwest`, which exists only in the native builds; the public router of web2 has no examples download route
- The example datasets are optional data that an installed app downloads from the GitHub release

> [!IMPORTANT]
> **Decision required.** How web2 offers example datasets:
>
> - the page downloads `examples.zip` of the release with `fetch` and the worker unpacks it into `/treetime/examples` (the release host must allow cross-origin requests)
> - the web2 deployment serves the small datasets of `data/` as static files next to the page
> - no examples in web2

## Validation

- The dataset catalog of web2 lists the example datasets, and a run of one of them computes
