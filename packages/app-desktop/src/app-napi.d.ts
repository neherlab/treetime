declare module "@neherlab/app-napi" {
  type EventCallback = (error: Error | null, eventJson: string) => void;

  type Command = (argsJson: string, onEvent: EventCallback) => Promise<string>;

  export function version(): string;
  export function datasets(): string;
  export function cancel(): void;
  export const ancestral: Command;
  export const clock: Command;
  export const timetree: Command;
  export const mugration: Command;
  export const optimize: Command;
  export const prune: Command;
}
