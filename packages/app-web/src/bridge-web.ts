import type { TreeTimeBridge, CommandResult, VersionInfo } from "@neherlab/app-contracts";

const API_BASE = "/api";

async function getJson<T>(path: string): Promise<T> {
  const response = await fetch(`${API_BASE}/${path}`);
  if (!response.ok) {
    throw new Error(`GET ${path}: ${response.status} ${response.statusText}`);
  }
  return (await response.json()) as T;
}

async function postCommand(command: string, args: unknown): Promise<CommandResult> {
  const response = await fetch(`${API_BASE}/${command}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(args),
  });
<<<<<<< HEAD

  if (!response.ok) {
    const text = await response.text();
    return { status: "error", error: text };
  }

  return response.json() as Promise<CommandResult>;
=======
  if (!response.ok) {
    return { status: "error", error: `${response.status} ${response.statusText}` };
  }
  return (await response.json()) as CommandResult;
>>>>>>> 33bee034 (feat: add end-to-end version info across all layers)
}

export function createWebBridge(): TreeTimeBridge {
  return {
    version: () => getJson<VersionInfo>("version"),
    ancestral: (args) => postCommand("ancestral", args),
    clock: (args) => postCommand("clock", args),
    timetree: (args) => postCommand("timetree", args),
    mugration: (args) => postCommand("mugration", args),
    optimize: (args) => postCommand("optimize", args),
    prune: (args) => postCommand("prune", args),
  };
}
