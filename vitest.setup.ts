import * as fc from "fast-check";

import { DETERMINISTIC_SEED } from "./test/seed";

fc.configureGlobal({ ...fc.readConfigureGlobal(), seed: DETERMINISTIC_SEED });

process.stderr.write(`[vitest] deterministic fast-check seed: ${DETERMINISTIC_SEED}\n`);
