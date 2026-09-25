export type * from "./generated/types.gen";

export * from "./generated/zod.gen";

export {
  BridgeError,
  bridgeErrorFromText,
  CancelledError,
  CommandError,
  RunEndedError,
  createBridge,
  parseRunEvent,
} from "./bridge";

export type {
  BridgeTransport,
  CheckConfigInput,
  CheckConfigResult,
  CladeInRunsResult,
  CommandOptions,
  DesktopRequestInput,
  FollowRunOptions,
  InputFactsResult,
  Parsed,
  RunConfigResult,
  RunComparisonResult,
  RunRecordResult,
  RunResultsResult,
  RunSummaryResult,
  TransportEventOptions,
  TreeTimeBridge,
} from "./bridge";

export { zPickedFiles, zPickFilesRequest, type LocalFiles, type PickFilesRequest } from "./files";

export { default as openApiDocument } from "../openapi.json";
