import type { TreeTimeBridge, CommandResult } from "@neherlab/app-contracts";

const API_BASE = "/api";

async function postCommand(command: string, args: unknown): Promise<CommandResult> {
  const response = await fetch(`${API_BASE}/${command}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(args),
  });

  if (!response.ok) {
    const text = await response.text();
    return { status: "error", error: text };
  }

  return response.json() as Promise<CommandResult>;
}

export function createWebBridge(): TreeTimeBridge {
  return {
    ancestral: (args) => postCommand("ancestral", args),
    clock: (args) => postCommand("clock", args),
    timetree: (args) => postCommand("timetree", args),
    mugration: (args) => postCommand("mugration", args),
    optimize: (args) => postCommand("optimize", args),
    prune: (args) => postCommand("prune", args),
  };
}
