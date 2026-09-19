export type CommandName = "ancestral" | "clock" | "mugration" | "optimize" | "prune" | "timetree";

export interface CommandInfo {
  name: CommandName;
  label: string;
  description: string;
}

export const COMMANDS: ReadonlyArray<CommandInfo> = [
  { name: "ancestral", label: "Ancestral", description: "Sequence reconstruction" },
  { name: "clock", label: "Clock", description: "Molecular clock test" },
  { name: "mugration", label: "Mugration", description: "Discrete trait reconstruction" },
  { name: "optimize", label: "Optimize", description: "Branch optimization" },
  { name: "prune", label: "Prune", description: "Tree pruning" },
  { name: "timetree", label: "Timetree", description: "Full timetree analysis" },
];

export const DEFAULT_COMMAND: CommandName = "ancestral";

export function isCommandName(value: string | undefined): value is CommandName {
  return COMMANDS.some((command) => command.name === value);
}
