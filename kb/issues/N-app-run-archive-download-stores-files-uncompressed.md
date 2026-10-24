# The downloaded run archive stores its files uncompressed

The web apps download the archive of a run's outputs through `GET /api/runs/{id}/archive`, which streams a zip whose entries are stored without compression. The server compresses the response on the wire, but the browser saves the decompressed zip, so the saved file is as large as the outputs. The `sc2/5100` timetree archive of web2 is 158 MB, most of it the reconstructed alignment. The desktop app writes the archive to a path with deflated entries instead.

## Evidence

- The streamed archive adds each file with `CompressionMethod::Stored` [packages/app-backend/src/runs/files.rs#L153-L156](../../packages/app-backend/src/runs/files.rs#L153-L156); the route streams it [packages/app-server/src/api/download.rs#L96](../../packages/app-server/src/api/download.rs#L96)
- `write_run_zip()`, which the desktop save uses, deflates each entry [packages/app-backend/src/runs/files.rs#L90-L92](../../packages/app-backend/src/runs/files.rs#L90-L92), [packages/app-backend/src/runs/manager.rs#L341](../../packages/app-backend/src/runs/manager.rs#L341)
- The server's compression layer skips only the content types it lists (xz, gzip, zstd, fonts, and the defaults of tower-http), so it compresses the stored zip on the wire [packages/app-server/src/web.rs#L138-L144](../../packages/app-server/src/web.rs#L138-L144); web2 has no wire, so its archive is not compressed at any step
- A route test fixes the stored entries [packages/app-server/src/__tests__/test_downloads.rs#L95](../../packages/app-server/src/__tests__/test_downloads.rs#L95)

> [!IMPORTANT]
> **Decision required.** The stored entries let the server stream the archive without compressing twice. Options: deflate the entries of the streamed archive (CPU per download, no compression on the wire needed for it), keep stored entries and accept saved archives as large as the outputs, or deflate only in web2, where nothing compresses the transfer.

## Validation

- The saved archive of a run matches the chosen compression, and the route test checks the chosen method
