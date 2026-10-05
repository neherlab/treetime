// oxlint-disable-next-line anti-slop/no-unknown-parameters -- a caught exception is untyped; this reads its message for display
export function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
