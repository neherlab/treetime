import * as fc from "fast-check";
import { inject } from "vitest";

fc.configureGlobal({ ...fc.readConfigureGlobal(), seed: inject("seed") });

declare module "vitest" {
  export interface ProvidedContext {
    seed: number;
  }
}
