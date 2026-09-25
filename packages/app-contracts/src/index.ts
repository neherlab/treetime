export type * from "./generated/types.gen";

export * from "./generated/zod.gen";

export { BridgeError, bridgeErrorFromText, CancelledError, CommandError, RunEndedError, createBridge } from "./bridge";

export type {
  BridgeTransport,
  DesktopRequestInput,
  InputFactsResult,
  Parsed,
  RunComparisonResult,
  RunRecordResult,
  RunResultsResult,
  RunSummaryResult,
  SettingKey,
  TransportEventOptions,
  TreeTimeBridge,
} from "./bridge";

export { zPickedFiles, zPickFilesRequest, type LocalFiles, type PickFilesRequest } from "./files";

export { default as openApiDocument } from "../openapi.json";
