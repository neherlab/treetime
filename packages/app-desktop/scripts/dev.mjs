import { createServer } from "vite";
import { spawn } from "child_process";
import { resolve, dirname } from "path";
import { fileURLToPath } from "url";

const pkgDir = resolve(dirname(fileURLToPath(import.meta.url)), "..");

const server = await createServer({
  configFile: resolve(pkgDir, "vite.config.ts"),
});
await server.listen();
const url = server.resolvedUrls?.local[0] ?? `http://localhost:${server.config.server.port}`;

const electronBin = resolve(pkgDir, "../../node_modules/.bin/electron");
const electron = spawn(electronBin, [".", "--no-sandbox"], {
  cwd: pkgDir,
  env: { ...process.env, VITE_DEV_SERVER_URL: url },
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
