import type { RunRecord } from "@neherlab/app-contracts";
import { configCheck } from "@neherlab/app-contracts/client";
import { keepPreviousData } from "@tanstack/react-query";
import { Fragment, useMemo, useState } from "react";
import RotateCcw from "~icons/lucide/rotate-ccw";

import { CodeLineView, keyedLines } from "../analysis/CodePanel";
import { defaultText } from "../analysis/SettingField";
import { useApi } from "../api/hooks";
import { CopyButton } from "../components/CopyButton";
import { OptionToggle } from "../components/OptionToggle";
import { Panel } from "../components/Panel";
import { COMMAND_SETTINGS, groupedSpecs } from "../settings/catalog";
import { settingValue } from "../settings/config";
import { zJsonObject } from "../settings/json";
import type { CodeFormat } from "../store/draftSchema";
import { Button } from "../ui/button";
import { cn } from "../ui/cn";
import { Label } from "../ui/label";
import { Switch } from "../ui/switch";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "../ui/table";
import { useRerun } from "./useRerun";

const CODE_FORMATS: ReadonlyArray<{ value: CodeFormat; label: string }> = [
  { value: "cli", label: "CLI" },
  { value: "yaml", label: "YAML" },
];

export function SettingsTab({ record }: { record: RunRecord }) {
  const settings = COMMAND_SETTINGS[record.command];
  const specs = settings.specs;
  const config = useMemo(() => zJsonObject.parse(record.config), [record.config]);
  const [changedOnly, setChangedOnly] = useState(true);
  const [format, setFormat] = useState<CodeFormat>("cli");
  const rerun = useRerun(record);
  const changedKeys = useMemo(() => new Set(record.changed_settings), [record.changed_settings]);

  const groups = useMemo(
    () =>
      groupedSpecs(settings, specs).flatMap(([group, members]) => {
        const rows = members.filter(
          (spec) => spec.role !== "output" && (!changedOnly || spec.role === "input" || changedKeys.has(spec.key)),
        );

        return rows.length === 0 ? [] : [{ group, rows }];
      }),
    [changedKeys, changedOnly, settings, specs],
  );

  const changedCount = changedKeys.size;

  const { data: check } = useApi(
    (context) =>
      configCheck({
        ...context,
        body: { command: record.command, text: JSON.stringify(config), input_facts: null },
      }),
    { placeholderData: keepPreviousData, staleTime: Infinity },
  );

  const code = check?.status === "valid" ? check.code : null;
  const lines = code === null ? [] : format === "cli" ? code.command_line : code.yaml;
  const text = code === null ? "" : format === "cli" ? code.command_line_text : code.yaml_text;

  return (
    <div className="grid gap-4 @5xl:grid-cols-[minmax(0,1fr)_32rem]">
      <Panel
        title="Settings used"
        hint={`${changedCount} ${changedCount === 1 ? "setting differs" : "settings differ"} from the defaults`}
        actions={
          <div className="flex items-center gap-2">
            <Switch id="settings-changed-only" checked={changedOnly} onCheckedChange={setChangedOnly} />
            <Label htmlFor="settings-changed-only">Changed only</Label>
          </div>
        }
      >
        <Table>
          <TableHeader>
            <TableRow>
              <TableHead>Setting</TableHead>
              <TableHead>Value</TableHead>
              <TableHead>Default</TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {groups.map(({ group, rows }) => (
              <Fragment key={group}>
                <TableRow className="bg-muted/50 hover:bg-muted/50">
                  <TableHead colSpan={3} scope="colgroup" className="h-7 text-xs">
                    {group}
                  </TableHead>
                </TableRow>
                {rows.map((spec) => {
                  const changed = changedKeys.has(spec.key);

                  return (
                    <TableRow key={spec.key} className={cn(changed && "bg-accent/60")}>
                      <TableCell className="whitespace-normal">
                        {changed && <span className="sr-only">Changed: </span>}
                        <span className={cn(changed && "font-bold")}>{spec.label}</span>{" "}
                        <code className="text-muted-foreground font-mono text-xs">{spec.flag}</code>
                      </TableCell>
                      <TableCell className="font-mono text-xs break-all whitespace-normal">
                        {defaultText(settingValue(config, spec))}
                      </TableCell>
                      <TableCell className="text-muted-foreground font-mono text-xs whitespace-normal">
                        {spec.role === "setting" ? defaultText(spec.default_value) : ""}
                      </TableCell>
                    </TableRow>
                  );
                })}
              </Fragment>
            ))}
          </TableBody>
        </Table>
      </Panel>
      <div className="grid content-start gap-4">
        <Panel
          title="Reproduce"
          actions={
            <>
              <OptionToggle label="Format" value={format} onChange={setFormat} options={CODE_FORMATS} />
              <CopyButton text={text} label={format === "cli" ? "Copy the command" : "Copy the config"} />
            </>
          }
        >
          <div className="grid gap-2 p-3.5">
            <pre className="bg-muted/50 overflow-x-auto rounded-md border px-3 py-2.5 font-mono text-xs leading-relaxed">
              {keyedLines(lines).map(({ key, line, index }) => (
                <CodeLineView
                  key={key}
                  line={line}
                  continued={format === "cli" && index < lines.length - 1}
                  indent={format === "cli" && index > 0}
                />
              ))}
            </pre>
            <p className="text-muted-foreground text-xs">
              TreeTime {record.treetime_version}. The command includes the outputs the app adds to every run and writes
              them to <code className="font-mono">out</code>.
            </p>
          </div>
        </Panel>
        <Button type="button" onClick={rerun}>
          <RotateCcw aria-hidden />
          Edit and run again
        </Button>
      </div>
    </div>
  );
}
