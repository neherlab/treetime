# The downloaded run archive stores its files uncompressed

web2 downloads the archive of a run's outputs through `GET /api/runs/{id}/archive`, which streams a zip whose entries are stored without compression. Nothing compresses the transfer from the back-end worker to the page, so the saved file is as large as the outputs. The `sc2/5100` timetree archive of web2 is 158 MB, most of it the reconstructed alignment. The desktop app writes the archive to a path with deflated entries instead.

## Evidence

- The streamed archive adds each file with `CompressionMethod::Stored` [packages/app-backend/src/runs/files.rs#L153-L156](../../packages/app-backend/src/runs/files.rs#L153-L156); the route streams it [packages/app-api/src/api/download.rs#L96](../../packages/app-api/src/api/download.rs#L96)
- `write_run_zip()`, which the desktop save uses, deflates each entry [packages/app-backend/src/runs/files.rs#L90-L92](../../packages/app-backend/src/runs/files.rs#L90-L92), [packages/app-backend/src/runs/manager.rs#L341](../../packages/app-backend/src/runs/manager.rs#L341)
- The router has no compression layer, and web2 passes the response through a message port, so its archive is not compressed at any step
- A route test fixes the stored entries [packages/app-api/src/__tests__/test_downloads.rs#L95](../../packages/app-api/src/__tests__/test_downloads.rs#L95)

> [!IMPORTANT]
> **Decision required.** Options: deflate the entries of the streamed archive (CPU per download, in the back-end worker that also answers the other requests of the tab), or keep stored entries and accept saved archives as large as the outputs.

## Validation

- The saved archive of a run matches the chosen compression, and the route test checks the chosen method
