import type { RunRecordResult } from "@neherlab/app-contracts";
import { Copy } from "lucide-react";
import { Fragment, useCallback, useMemo, useState } from "react";

import { CodeLineView, keyedLines } from "../analysis/CodePanel";
import { defaultText } from "../analysis/SettingField";
import { useConfigCheck } from "../queries";
import { COMMAND_SETTINGS, groupedSpecs } from "../settings/catalog";
import { outputFreeConfig, settingValue } from "../settings/config";
import { zJsonObject } from "../settings/json";
import type { CodeFormat } from "../store/draftSchema";
import { Button, Segmented, Switch, cn } from "../ui";
import { Panel } from "./Panel";
import { useCopy } from "./useCopy";
import { useRerun } from "./useRerun";

const CODE_FORMATS: ReadonlyArray<{ value: CodeFormat; label: string }> = [
  { value: "cli", label: "CLI" },
  { value: "yaml", label: "YAML" },
];

export function SettingsTab({ record }: { record: RunRecordResult }) {
  const settings = COMMAND_SETTINGS[record.command];
  const specs = settings.specs;
  const config = useMemo(() => zJsonObject.parse(record.config), [record.config]);
  const [changedOnly, setChangedOnly] = useState(true);
  const [format, setFormat] = useState<CodeFormat>("cli");
  const copy = useCopy();
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
  const reproducible = useMemo(() => outputFreeConfig(specs, config), [config, specs]);
  const { data: check } = useConfigCheck(record.command, reproducible, null);
  const code = check?.status === "valid" ? check.code : null;
  const lines = code === null ? [] : format === "cli" ? code.command_line : code.yaml;
  const text = code === null ? "" : format === "cli" ? code.command_line_text : code.yaml_text;

  const onCopy = useCallback(
    () => copy(text, format === "cli" ? "Command copied to the clipboard" : "Config copied to the clipboard"),
    [copy, format, text],
  );

  return (
    <div className="grid gap-3.5 xl:grid-cols-[minmax(0,1fr)_32rem]">
      <Panel
        title="Settings used"
        hint={`${changedCount} ${changedCount === 1 ? "setting differs" : "settings differ"} from the defaults`}
        actions={
          <span className="text-ink-muted inline-flex items-center gap-2 text-sm">
            <Switch checked={changedOnly} onCheckedChange={setChangedOnly} label="Changed only" />
            <span aria-hidden>Changed only</span>
          </span>
        }
      >
        <table className="w-full border-collapse text-left">
          <thead>
            <tr className="text-ink-faint text-xs">
              <th className="px-3.5 py-1.5 font-normal">Setting</th>
              <th className="px-3.5 py-1.5 font-normal">Value</th>
              <th className="px-3.5 py-1.5 font-normal">Default</th>
            </tr>
          </thead>
          <tbody>
            {groups.map(({ group, rows }) => (
              <Fragment key={group}>
                <tr className="bg-surface-2">
                  <th colSpan={3} className="text-ink-muted px-3.5 py-1 text-xs font-bold">
                    {group}
                  </th>
                </tr>
                {rows.map((spec) => {
                  const changed = changedKeys.has(spec.key);

                  return (
                    <tr key={spec.key} className={cn("border-line border-t", changed && "bg-accent-subtle")}>
                      <td className="px-3.5 py-1.5">
                        {changed && <span className="sr-only">Changed: </span>}
                        <span className={cn(changed && "font-bold")}>{spec.label}</span>{" "}
                        <code className="text-ink-faint font-mono text-xs">{spec.flag}</code>
                      </td>
                      <td className="px-3.5 py-1.5 font-mono text-xs break-all">
                        {defaultText(settingValue(config, spec))}
                      </td>
                      <td className="text-ink-faint px-3.5 py-1.5 font-mono text-xs">
                        {spec.role === "setting" ? defaultText(spec.default_value) : ""}
                      </td>
                    </tr>
                  );
                })}
              </Fragment>
            ))}
          </tbody>
        </table>
      </Panel>
      <div className="grid content-start gap-3.5">
        <Panel
          title="Reproduce"
          actions={
            <>
              <Segmented label="Format" value={format} onChange={setFormat} options={CODE_FORMATS} />
              <Button type="button" variant="ghost" size="icon" aria-label="Copy" onClick={onCopy}>
                <Copy size={14} aria-hidden />
              </Button>
            </>
          }
        >
          <div className="grid gap-2 px-3.5 py-3">
            <pre className="border-line bg-surface-2 m-0 max-h-96 overflow-auto rounded-md border px-3 py-2.5 font-mono text-xs leading-relaxed">
              {keyedLines(lines).map(({ key, line, index }) => (
                <CodeLineView
                  key={key}
                  line={line}
                  continued={format === "cli" && index < lines.length - 1}
                  indent={format === "cli" && index > 0}
                />
              ))}
            </pre>
            <p className="text-ink-faint m-0 text-xs">
              TreeTime {record.treetime_version}. The command includes the outputs the app adds to every run and writes
              them to <code className="font-mono">out</code>.
            </p>
          </div>
        </Panel>
        <Button type="button" onClick={rerun}>
          Edit and run again
        </Button>
      </div>
    </div>
  );
}
