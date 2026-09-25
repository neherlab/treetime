const CRASH_WINDOW_MS = 60_000;

const MAX_CRASHES_IN_WINDOW = 5;

export function shouldRestart(exitTimes: readonly number[], now: number): boolean {
  return exitTimes.filter((time) => now - time < CRASH_WINDOW_MS).length < MAX_CRASHES_IN_WINDOW;
}
