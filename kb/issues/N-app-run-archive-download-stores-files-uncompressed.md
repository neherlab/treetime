# The downloaded run archive stores its files uncompressed

The web apps download the archive of a run's outputs through `GET /api/runs/{id}/archive`, which streams a zip whose entries are stored without compression. The server compresses the response on the wire, but the browser saves the decompressed zip, so the saved file is as large as the outputs. The `sc2/5100` timetree archive of web2 is 158 MB, most of it the reconstructed alignment. The desktop app writes the archive to a path with deflated entries instead.

## Evidence

- The streamed archive adds each file with `CompressionMethod::Stored` [packages/app-backend/src/runs/files.rs#L151-L156](../../packages/app-backend/src/runs/files.rs#L151-L156); the route streams it [packages/app-server/src/api/download.rs#L97](../../packages/app-server/src/api/download.rs#L97)
- `write_run_zip()`, which the desktop save uses, deflates each entry [packages/app-backend/src/runs/files.rs#L90-L92](../../packages/app-backend/src/runs/files.rs#L90-L92), [packages/app-backend/src/runs/manager.rs#L348](../../packages/app-backend/src/runs/manager.rs#L348)
- The server's compression layer skips only the content types it lists (xz, gzip, zstd, fonts, and the defaults of tower-http), so it compresses the stored zip on the wire [packages/app-server/src/web.rs#L138-L144](../../packages/app-server/src/web.rs#L138-L144); web2 has no wire, so its archive is not compressed at any step
- A route test fixes the stored entries [packages/app-server/src/__tests__/test_downloads.rs#L98](../../packages/app-server/src/__tests__/test_downloads.rs#L98)

> [!IMPORTANT]
> **Decision required.** The stored entries let the server stream the archive without compressing twice. Options: deflate the entries of the streamed archive (CPU per download, no compression on the wire needed for it), keep stored entries and accept saved archives as large as the outputs, or deflate only in web2, where nothing compresses the transfer.

## Validation

- The saved archive of a run matches the chosen compression, and the route test checks the chosen method

## Proposed implementation plan

Every cited location is current (`fn stream_run_zip()` sets `CompressionMethod::Stored` at [packages/app-backend/src/runs/files.rs#L155](../../packages/app-backend/src/runs/files.rs#L155)). The 158 MB size of the `sc2/5100` archive was not measured again. The plan stops at the open decision above; one short plan per option:

- **Deflate the streamed archive**: set `CompressionMethod::Deflated` in `fn stream_run_zip()` as `fn write_run_zip()` does, and add `application/zip` to `fn compress_predicate()` in `packages/app-server/src/web.rs` as a content type that the wire layer skips, so the server does not compress twice. Update the expected method in `packages/app-server/src/__tests__/test_downloads.rs#L98`
- **Keep stored entries**: no code change; delete this issue after the decision
- **Deflate only in web2**: pass the compression method as a parameter of `fn stream_run_zip()` from its callers (the server route keeps `Stored`, the WebAssembly back end passes `Deflated`), and add a test of the method in the wasm back end tests

For a deflate option, a backend test in `packages/app-backend/src/runs/__tests__/` streams a run with a repetitive output file and asserts that every entry's method is `Deflated` and that the unpacked bytes equal the file. Oracle: the zip format (entry method field) and the input file.

### Plan validation

- `just lint-rs`, `just test-rs`, and for the web2 option `just test-wasm` pass
