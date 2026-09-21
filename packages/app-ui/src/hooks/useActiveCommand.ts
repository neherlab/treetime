import { useParams } from "@tanstack/react-router";

import { DEFAULT_COMMAND, isCommandName, type CommandName } from "../types";

export function useActiveCommand(): CommandName {
  const params = useParams({ strict: false });

  return isCommandName(params.command) ? params.command : DEFAULT_COMMAND;
}
