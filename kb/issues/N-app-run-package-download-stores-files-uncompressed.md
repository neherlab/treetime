# The downloaded run package stores its files uncompressed

The page downloads the run package of a run through `GET /api/runs/{id}/package`, which streams a zip whose entries are stored without compression. Nothing compresses the transfer from the back-end worker to the page, so the saved package is as large as its inputs and outputs together. The `sc2/5100` timetree outputs of web2 alone are 158 MB, most of it the reconstructed alignment, and the largest example package (ancestral on `mpox/clade-ii/2000`, default outputs) is about 1.18 GB. The desktop app writes the package to a path with deflated entries instead.

## Evidence

- The streamed package adds each entry with `CompressionMethod::Stored` (`RunPackageChunks` in [packages/app-backend/src/runs/package.rs#L183-L186](../../packages/app-backend/src/runs/package.rs#L183-L186)); the route streams it (`stream_run_package()`, [packages/app-api/src/api/download.rs#L65](../../packages/app-api/src/api/download.rs#L65))
- `RunPackage::write()`, which the desktop save uses, deflates each entry [packages/app-backend/src/runs/package.rs#L118-L120](../../packages/app-backend/src/runs/package.rs#L118-L120), [packages/app-backend/src/runs/manager.rs#L486](../../packages/app-backend/src/runs/manager.rs#L486)
- The router has no compression layer, and web2 passes the response through a message port, so its package is not compressed at any step
- A route test fixes the stored entries [packages/app-api/src/__tests__/test_downloads.rs#L88](../../packages/app-api/src/__tests__/test_downloads.rs#L88)

> [!IMPORTANT]
> **Decision required.** Options: deflate the entries of the streamed package (CPU per download, in the back-end worker that also answers the other requests of the tab), or keep stored entries and accept packages as large as the run files.

## Validation

- The saved package of a run matches the chosen compression, and the route test checks the chosen method
