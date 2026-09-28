import type { AnyAction } from "redux";
import type { ThunkAction } from "redux-thunk";

export interface AuspiceState {
  metadata: AuspiceMetadataState;
  tree: AuspiceTreeState;
  controls: AuspiceControlsState;
}

export type AuspiceThunk = ThunkAction<void, AuspiceState, undefined, AnyAction>;

export interface AuspiceMetadataState {
  loaded: boolean;
}

export interface AuspiceTreeState {
  loaded: boolean;
  nodes: AuspiceNode[] | null;
  idxOfInViewRootNode: number;
  visibility: number[] | null;
}

interface AuspiceNode {
  name: string;
  arrayIdx: number;
  hasChildren: boolean;
}

export interface AuspiceControlsState {
  colorBy: string;
  distanceMeasure: string;
  panelsToDisplay: readonly string[];
  selectedNode: AuspiceSelectedNode | null;
  filters: Readonly<Record<string, readonly AuspiceFilterValue[]>>;
  performanceFlags: ReadonlyMap<string, boolean>;
  scatterVariables: AuspiceScatterVariables;
}

interface AuspiceFilterValue {
  value: string;
  active: boolean;
}

interface AuspiceScatterVariables {
  showBranches?: boolean;
  showRegression?: boolean;
  x?: string;
  xLabel?: string;
  xContinuous?: boolean;
  xDomain?: readonly (string | number)[];
  xTemporal?: boolean;
  y?: string;
  yLabel?: string;
  yContinuous?: boolean;
  yDomain?: readonly (string | number)[];
  yTemporal?: boolean;
}

interface AuspiceSelectedNode {
  name: string;
  idx: number;
  isBranch: boolean;
}

export interface AuspicePublication {
  author: string;
  title: string;
  journal: string;
  year: string;
  href: string;
}
