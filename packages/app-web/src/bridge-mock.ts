import type { TreeTimeBridge, CommandResult, VersionInfo, ProgressEvent } from "@neherlab/app-contracts";

function delay(ms: number): Promise<void> {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

const PROGRESS_PHASES = [
  { stage: "Parsing input", fraction: 0.1 },
  { stage: "Building tree", fraction: 0.25 },
  { stage: "Inferring GTR model", fraction: 0.4 },
  { stage: "Optimizing branch lengths", fraction: 0.6 },
  { stage: "Iteration 2 of 3", fraction: 0.75 },
  { stage: "Ancestral reconstruction", fraction: 0.9 },
  { stage: "Writing output", fraction: 1.0 },
];

async function simulateRun(onProgress?: (event: ProgressEvent) => void): Promise<CommandResult> {
  for (const phase of PROGRESS_PHASES) {
    await delay(400);
    onProgress?.({ stage: phase.stage, fraction: phase.fraction, message: phase.stage });
  }
  return { status: "ok" };
}

export function createMockBridge(): TreeTimeBridge {
  return {
    version: async (): Promise<VersionInfo> => ({ version: "0.1.0-mockup" }),
    ancestral: (_args, onProgress) => simulateRun(onProgress),
    clock: (_args, onProgress) => simulateRun(onProgress),
    timetree: (_args, onProgress) => simulateRun(onProgress),
    mugration: (_args, onProgress) => simulateRun(onProgress),
    optimize: (_args, onProgress) => simulateRun(onProgress),
    prune: (_args, onProgress) => simulateRun(onProgress),
  };
}
