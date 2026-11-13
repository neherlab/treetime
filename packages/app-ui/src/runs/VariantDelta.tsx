import type { RunRecord } from "@neherlab/app-contracts";
import { runsCompare } from "@neherlab/app-contracts/client";
import { Link } from "@tanstack/react-router";

import { useApi } from "../api/hooks";

export function VariantDelta({ record }: { record: RunRecord }) {
  const source = record.source_run_id;

  const { data: comparison } = useApi(
    (context) => runsCompare({ ...context, path: { id: source ?? "", other: record.id } }),
    { enabled: source !== undefined, staleTime: Infinity },
  );

  const summary = comparison?.delta?.summary;

  if (source === undefined || summary === undefined) {
    return null;
  }

  return (
    <p className="flex flex-wrap items-baseline gap-x-3 gap-y-1 text-sm">
      <span className="tabular-nums">{summary} against the source</span>
      <Link
        to="/compare/$a/$b"
        params={{ a: source, b: record.id }}
        className="text-primary font-bold underline-offset-4 hover:underline"
      >
        Full comparison
      </Link>
    </p>
  );
}
