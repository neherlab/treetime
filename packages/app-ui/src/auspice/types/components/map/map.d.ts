import type { ComponentType } from "react";

declare const Map: ComponentType<{
  width: number;
  height: number;
  justGotNewDatasetRenderNewMap: boolean;
  legend: { swatches: boolean; demes: boolean };
}>;

export default Map;
