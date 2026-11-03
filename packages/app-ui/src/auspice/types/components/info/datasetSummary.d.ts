import type { AuspiceNode, AuspiceTranslate } from "../../../state";

export declare function datasetSummary(options: {
  nodes: readonly AuspiceNode[];
  visibility: readonly number[];
  mainTreeNumTips: number;
  temporalTree: boolean;
  t: AuspiceTranslate;
}): string;
