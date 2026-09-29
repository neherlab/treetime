import { fileURLToPath } from "node:url";

export const publicDir = fileURLToPath(new URL("../public", import.meta.url));
