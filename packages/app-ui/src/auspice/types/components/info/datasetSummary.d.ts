import type { AuspiceNode, AuspiceTranslate } from "../../../state";

export declare function datasetSummary(options: {
  nodes: readonly AuspiceNode[];
  visibility: readonly number[];
  mainTreeNumTips: number;
  branchLengthsToDisplay: string;
  t: AuspiceTranslate;
}): string;
