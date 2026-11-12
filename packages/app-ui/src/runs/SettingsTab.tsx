import type { RunRecord, SparseConfig } from "@neherlab/app-contracts";
import { Label } from "@neherlab/design/ui/label";
import { Switch } from "@neherlab/design/ui/switch";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/design/ui/table";
import { cn } from "cn";
import { Fragment, useMemo, useState } from "react";

import { Panel } from "../components/Panel";
import { useLocale } from "../i18n/locale";
import { commandSettings } from "../settings/catalog";
import { documentedDefault, effectiveValue } from "../settings/config";
import { groupedSpecs } from "../settings/groups";
import { defaultText } from "../settings/help";
import { shownSettings } from "./settingRows";

export function SettingsTab({ record }: { record: RunRecord }) {
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

  const changedCount = changedKeys.size;

  return (
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
                      {defaultText(effectiveValue(config, spec), locale)}
                    </TableCell>
                    <TableCell className="text-muted-foreground font-mono text-xs whitespace-normal">
                      {spec.role === "setting" ? defaultText(documentedDefault(spec), locale) : ""}
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
