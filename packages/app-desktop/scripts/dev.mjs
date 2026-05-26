import { createServer } from "vite";
import { spawn } from "child_process";
import { resolve, dirname } from "path";
import { fileURLToPath } from "url";
import { existsSync } from "fs";

const pkgDir = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const napiDir = resolve(pkgDir, "../app-napi");

const server = await createServer({
  configFile: resolve(pkgDir, "vite.config.ts"),
});
await server.listen();
const url = server.resolvedUrls?.local[0] ?? `http://localhost:${server.config.server.port}`;

const napiNode = resolve(napiDir, "app-napi.linux-x64-gnu.node");
const ldPreload = [process.env.LD_PRELOAD, existsSync(napiNode) ? napiNode : ""]
  .filter(Boolean)
  .join(":");

const projectRoot = resolve(pkgDir, "../..");
const electronBin = resolve(projectRoot, "node_modules/.bin/electron");
const electron = spawn(electronBin, [".", "--no-sandbox"], {
  cwd: pkgDir,
  env: {
    ...process.env,
    VITE_DEV_SERVER_URL: url,
    TREETIME_PROJECT_ROOT: projectRoot,
    ...(ldPreload ? { LD_PRELOAD: ldPreload } : {}),
  },
  stdio: "inherit",
});

electron.on("close", () => {
  server.close();
  process.exit();
});

process.on("SIGINT", () => {
  electron.kill();
  server.close();
  process.exit();
});
