---
title: Developer guide
description: How to set up, build, test, and release TreeTime, and how the repository is organized.
---

This guide describes how to set up a development environment, build and test TreeTime, and maintain the project. It assumes basic familiarity with TreeTime from a user perspective.

The `justfile` is the entry point for every routine task. `just` lists the recipes by group, and each recipe has a one-line description.

## Project structure

```
packages/
  app-backend/       Run records, results, settings, and pre-run checks of the server, the addon, and the web2 worker
  app-cli/           CLI binary (treetime)
  app-commands/      Command configs and runners shared by the CLI, the server, and the addon
  app-contracts/     TypeScript types and schemas of the bridge between the apps and the Rust code
  app-datasets/      Bundled example datasets
  app-desktop/       Electron shell of the desktop app
  app-docs/          Docs site code (Astro with headless fumadocs-core)
  app-napi/          Node addon of the desktop app (Rust, napi-rs)
  app-output/        Output writers shared by the CLI, the server, and the addon
  app-server/        HTTP API server of the web app (treetime-server)
  app-ui/            React components shared by the desktop and web apps
  app-wasm/          WebAssembly back end of the web2 app (Rust, Emscripten, wasm-bindgen)
  app-web/           Web app (Vite)
  app-web2/          Web app on the WebAssembly back end in a web worker (Vite)
  legacy/            TreeTime v0, the Python reference implementation
  schemas/           Generated JSON schemas of the CLI input
  treetime/          Core library
  treetime-*/        Supporting crates of the core library
  util-*/            File format libraries
dev/                 Development scripts and the container setup
docs/                Pages of the docs site, one folder per language
test_scripts/        Python research scripts and notebooks
kb/                  Knowledge base
data/                Example datasets
```
