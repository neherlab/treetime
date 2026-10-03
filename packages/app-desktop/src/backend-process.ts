const CRASH_WINDOW_MS = 60_000;

const MAX_CRASHES_IN_WINDOW = 5;

export function shouldRestart(exitTimes: readonly number[], now: number): boolean {
  return exitTimes.filter((time) => now - time < CRASH_WINDOW_MS).length < MAX_CRASHES_IN_WINDOW;
}

export function stopReason(code: number, requested: boolean, restart: boolean): string {
  if (requested) {
    return "the back end restarts to open the new runs folder; the request was not answered";
  }

  return restart
    ? `the back end stopped with exit code ${code} and restarts; the request was not answered`
    : `the back end stopped with exit code ${code} too often and does not restart; restart TreeTime`;
}
