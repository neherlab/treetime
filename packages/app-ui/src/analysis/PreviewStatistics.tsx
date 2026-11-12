import type { RegressionStatistic } from "@neherlab/app-contracts";
import { HoverCard, HoverCardContent, HoverCardTrigger } from "@neherlab/design/ui/hover-card";

import { statisticText } from "./draftPreview";

export function PreviewStatistics({ statistics }: { statistics: readonly RegressionStatistic[] }) {
  if (statistics.length === 0) {
    return null;
  }

  return (
    <dl className="flex flex-wrap gap-x-6 gap-y-2 px-2 pb-2 text-sm">
      {statistics.map((statistic) => (
        <div key={statistic.statistic} className="grid gap-0.5">
          <dt className="text-muted-foreground text-xs">
            <HoverCard>
              <HoverCardTrigger
                render={
                  <button
                    type="button"
                    aria-label={statistic.label}
                    className="cursor-help underline decoration-dotted underline-offset-4"
                  />
                }
              >
                {statistic.label}
              </HoverCardTrigger>
              <HoverCardContent className="w-72 text-sm">{statistic.definition}</HoverCardContent>
            </HoverCard>
          </dt>
          <dd className="font-bold tabular-nums">{statisticText(statistic)}</dd>
        </div>
      ))}
    </dl>
  );
}
