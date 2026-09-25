import type { RunRecordResult } from "@neherlab/app-contracts";
import { Link } from "@tanstack/react-router";

import { defaultText } from "../analysis/SettingField";
import { useRunRecords } from "../queries";
import { COMMAND_SETTINGS } from "../settings/catalog";
import { COMMAND_INFO } from "../settings/commands";
import { settingValue } from "../settings/config";
import { sameJson, zJsonObject } from "../settings/json";
import { settingLabel } from "../settings/labels";

export function ComparePage({ first, second }: { first: string; second: string }) {
  const [left, right] = useRunRecords([first, second]);

  if (left?.data === undefined || right?.data === undefined) {
    const error = left?.error ?? right?.error ?? null;

    return (
      <p className="text-ink-muted p-10 text-center">
        {error === null ? "Loading the runs..." : `The runs cannot be loaded: ${error.message}`}
      </p>
    );
  }

  return <Comparison left={left.data} right={right.data} />;
}

function Comparison({ left, right }: { left: RunRecordResult; right: RunRecordResult }) {
  const leftConfig = zJsonObject.parse(left.config);
  const rightConfig = zJsonObject.parse(right.config);
  const specs = left.command === right.command ? COMMAND_SETTINGS[left.command].specs : [];

  const differing = specs.filter(
    (spec) => spec.pathRole !== "output" && !sameJson(settingValue(leftConfig, spec), settingValue(rightConfig, spec)),
  );

  return (
    <div className="mx-auto max-w-[92.5rem] px-5 pt-4 pb-16">
      <h1 className="mb-3.5 text-2xl font-bold">Compare runs</h1>
      <div className="border-line bg-surface-1 rounded-lg border">
        <table className="w-full border-collapse">
          <thead>
            <tr className="text-left">
              <th className="text-ink-faint px-3.5 py-2 text-xs">Setting</th>
              {[left, right].map((record) => (
                <th key={record.id} className="px-3.5 py-2">
                  <Link to="/runs/$id/results" params={{ id: record.id }} className="text-accent">
                    {record.title}
                  </Link>
                  <span className="text-ink-faint block text-xs font-normal">{COMMAND_INFO[record.command].label}</span>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {left.command !== right.command && (
              <tr className="border-line border-t">
                <td colSpan={3} className="text-ink-muted px-3.5 py-2">
                  The runs use different analyses, so their settings are not compared.
                </td>
              </tr>
            )}
            {left.command === right.command && differing.length === 0 && (
              <tr className="border-line border-t">
                <td colSpan={3} className="text-ink-muted px-3.5 py-2">
                  The runs have the same settings and input paths.
                </td>
              </tr>
            )}
            {differing.map((spec) => (
              <tr key={spec.key} className="border-line border-t">
                <td className="px-3.5 py-1.5">
                  {settingLabel(spec.key)} <code className="text-ink-faint font-mono text-xs">{spec.flag}</code>
                </td>
                <td className="px-3.5 py-1.5 font-mono text-xs">{defaultText(settingValue(leftConfig, spec))}</td>
                <td className="px-3.5 py-1.5 font-mono text-xs">{defaultText(settingValue(rightConfig, spec))}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
