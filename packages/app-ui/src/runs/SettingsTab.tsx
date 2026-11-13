import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";
import { Button } from "@neherlab/design/ui/button";
import { Label } from "@neherlab/design/ui/label";
import { Switch } from "@neherlab/design/ui/switch";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";
import { cn } from "cn";
import { Fragment, useCallback, useMemo, useState } from "react";
import Download from "~icons/lucide/download";

import { Panel } from "../components/Panel";
import { useLocale } from "../i18n/locale";
import { commandSettings } from "../settings/catalog";
import { documentedDefault, effectiveValue } from "../settings/config";
import { groupedSpecs } from "../settings/groups";
import { defaultText, SettingHelpPopover } from "../settings/help";
import { useRunCode } from "./RunCommandLine";
import { saveText } from "./saveRunOutput";
import { changedSettingCount, shownSettings } from "./settingRows";

export function SettingsTab({ record, configName }: { record: RunRecord; configName: string | undefined }) {
  const locale = useLocale();
  const settings = commandSettings(record.command);
  const specs = settings.settings;
  const config: SparseConfig = record.config;
  const [changedOnly, setChangedOnly] = useState(true);
  const changedKeys = useMemo(() => new Set(record.changed_settings), [record.changed_settings]);

  const groups = useMemo(
    () =>
      groupedSpecs(settings, specs).flatMap(([group, members]) => {
        const rows = shownSettings(members, config, changedKeys, changedOnly);

        return rows.length === 0 ? [] : [{ group, rows }];
      }),
    [changedKeys, changedOnly, config, settings, specs],
  );

  const changedCount = changedSettingCount(specs, changedKeys);
  const code = useRunCode(record, true);
  const yaml = code.state === "ready" ? code.code.yaml_text : undefined;
  const commandLines = code.state === "ready" ? code.commandLines : undefined;

  const onExport = useCallback(() => {
    if (yaml !== undefined && configName !== undefined) {
      saveText(yaml, configName, "application/yaml");
    }
  }, [configName, yaml]);

  return (
    <Panel
      title="Settings used"
      hint={`${changedCount} ${changedCount === 1 ? "setting differs" : "settings differ"} from the defaults`}
      actions={
        <div className="flex items-center gap-2">
          <Switch id="settings-changed-only" checked={changedOnly} onCheckedChange={setChangedOnly} />
          <Label htmlFor="settings-changed-only">Changed only</Label>
          <Button
            type="button"
            variant="outline"
            size="sm"
            onClick={onExport}
            disabled={yaml === undefined || configName === undefined}
          >
            <Download aria-hidden />
            Config (.yaml)
          </Button>
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
                      <SettingHelpPopover
                        spec={spec}
                        commandLine={commandLines?.find((line) => line.key === spec.key)}
                      />
                    </TableCell>
                    <TableCell className="font-mono text-xs break-all whitespace-normal">
                      {defaultText(effectiveValue(config, spec), locale, spec.options)}
                    </TableCell>
                    <TableCell className="text-muted-foreground font-mono text-xs whitespace-normal">
                      {spec.role === "setting" ? defaultText(documentedDefault(spec), locale, spec.options) : ""}
                    </TableCell>
                  </TableRow>
                );
              })}
            </Fragment>
          ))}
        </TableBody>
      </Table>
    </Panel>
  );
}
