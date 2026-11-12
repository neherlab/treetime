import type { RunRecord } from "@neherlab/app-contracts";
import { runConfig, runsSettings } from "@neherlab/app-contracts/client";
import { Button } from "@neherlab/design/ui/button";
import { useCallback, useMemo, useState } from "react";
import SquareTerminal from "~icons/lucide/square-terminal";

import { type CommandLineCode, CommandLineDialog } from "../analysis/CommandLineDialog";
import { useApi } from "../api/hooks";
import { rerunDraft } from "../settings/rerun";

export function RunCommandLine({ record }: { record: RunRecord }) {
  const [open, setOpen] = useState(false);
  const show = useCallback(() => setOpen(true), []);
  const code = useRunCode(record, open);

  const seeded = record.seed === undefined ? "" : `, seed ${record.seed}`;

  return (
    <>
      <Button type="button" variant="outline" size="sm" onClick={show}>
        <SquareTerminal aria-hidden />
        Command line
      </Button>
      <CommandLineDialog
        open={open}
        onOpenChange={setOpen}
        description={`${record.title}: TreeTime ${record.treetime_version}${seeded}`}
        code={code}
      />
    </>
  );
}

export function useRunCode(record: RunRecord, enabled: boolean): CommandLineCode {
  const { data: settings, error: settingsError } = useApi(
    (context) => runsSettings({ ...context, path: { id: record.id } }),
    { enabled, staleTime: Infinity },
  );

  const config = useMemo(
    () => (settings === undefined ? undefined : rerunDraft(record, settings).config),
    [record, settings],
  );

  const { data: resolved, error: resolveError } = useApi(
    (context) => runConfig({ ...context, body: { command: record.command, config: config ?? {} } }),
    { enabled: enabled && config !== undefined, staleTime: Infinity },
  );

  const error = settingsError ?? resolveError;

  if (error !== null) {
    return { state: "invalid", message: error.message };
  }

  if (resolved === undefined) {
    return { state: "loading" };
  }

  return resolved.status === "valid"
    ? { state: "ready", code: resolved.code }
    : { state: "invalid", message: [resolved.message, ...resolved.causes].join(": ") };
}
