import type { AnyAction } from "redux";
import type { ThunkAction } from "redux-thunk";

export interface AuspiceState {
  metadata: AuspiceMetadataState;
  tree: AuspiceTreeState;
  controls: AuspiceControlsState;
  entropy: AuspiceEntropyState;
  notifications: AuspiceNotificationsState;
  treeToo: AuspiceTreeTooState;
}

export type AuspiceThunk = ThunkAction<void, AuspiceState, undefined, AnyAction>;

export type AuspiceTranslate = (key: string, values?: Readonly<Record<string, string | number>>) => string;

export interface AuspiceMetadataState {
  loaded: boolean;
  title?: string;
  updated?: string;
  mainTreeNumTips: number;
  sharing: AuspiceSharing;
  colorings: Readonly<Record<string, AuspiceColoring>>;
}

interface AuspiceColoring {
  title?: string;
  type: string;
}

export interface AuspiceSharing {
  dataset_json?: boolean;
  metadata_tsv?: boolean;
  trees?: boolean;
  entropy?: boolean;
  screenshot?: boolean;
}

export interface AuspiceTreeState {
  loaded: boolean;
  nodes: AuspiceNode[] | null;
  idxOfInViewRootNode: number;
  visibility: number[] | null;
  nodeColors: readonly string[];
}

export interface AuspiceTreeTooState {
  loaded: boolean;
  nodes: AuspiceNode[] | null;
  nodeColors: readonly string[];
  tangleTipLookup?: readonly (readonly number[])[];
}

export interface AuspiceNode {
  name: string;
  arrayIdx: number;
  hasChildren: boolean;
  fullTipCount: number;
  parent: AuspiceNode;
}

export interface AuspiceControlsState {
  colorBy: string;
  distanceMeasure: string;
  branchLengthsToDisplay: string;
  panelLayout: AuspicePanelLayout;
  canTogglePanelLayout: boolean;
  sidebarOpen: boolean;
  showTreeToo: string | undefined;
  explodeAttr: string | undefined;
  panelsAvailable: readonly string[];
  panelsToDisplay: readonly string[];
  selectedNode: AuspiceSelectedNode | null;
  filters: Readonly<Partial<Record<string | symbol, readonly AuspiceFilterValue[]>>>;
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

export type AuspicePanelLayout = "full" | "grid";

export interface AuspiceEntropyState {
  loaded: boolean;
}

export interface AuspiceNotificationsState {
  stack: readonly AuspiceNotification[];
}

export interface AuspiceNotification {
  id: number;
  notificationType: AuspiceNotificationType;
  message: string;
  details: string;
}

export type AuspiceNotificationType = "info" | "success" | "warning" | "error";

export interface AuspicePublication {
  author: string;
  title: string;
  journal: string;
  year: string;
  href: string;
}
