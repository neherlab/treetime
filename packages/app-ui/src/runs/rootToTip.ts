import type { RootToTip } from "@neherlab/app-contracts";
import { useMemo } from "react";

import { formatRate } from "../format";
import type { RttLine, RttPoint } from "./RootToTipPlot";

export function useRootToTip(regression: RootToTip | undefined): {
  points: RttPoint[];
  line: RttLine | undefined;
} {
  return useMemo(
    () =>
      regression === undefined
        ? { points: [], line: undefined }
        : { points: rttPoints(regression), line: rttLine(regression) },
    [regression],
  );
}

function rttPoints(regression: RootToTip): RttPoint[] {
  return regression.points.flatMap((point) =>
    point.date === undefined
      ? []
      : [
          {
            name: point.name,
            date: point.date.year,
            dateText: point.date.date,
            div: point.div,
            excluded: point.outlier,
            inferred: point.date_source === "inferred",
          },
        ],
  );
}

function rttLine(regression: RootToTip): RttLine | undefined {
  const line = regression.line;

  return line === undefined
    ? undefined
    : {
        slope: line.rate,
        intercept: line.intercept,
        label: `TreeTime clock model: rate ${formatRate(line.rate)} /site/yr`,
      };
}
