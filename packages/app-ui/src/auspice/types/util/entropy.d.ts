import type { AuspiceEntropyBar, AuspiceNode, AuspiceSelectedCds } from "../../state";

export declare function calcEntropyInView(
  nodes: readonly AuspiceNode[],
  visibility: readonly number[],
  selectedCds: AuspiceSelectedCds,
  showCounts: boolean,
): [AuspiceEntropyBar[], number];
