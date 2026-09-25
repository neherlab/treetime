import type { ColorScales } from "../results/colors";
import { getAt, isJsonObject, isString, setAt, type JsonObject } from "../settings/json";

const COLORINGS_PATH = ["meta", "colorings"];

const SCALE_PATH = ["scale"];

export function withColorScales(document: JsonObject, scales: ColorScales): JsonObject {
  const colorings = getAt(document, COLORINGS_PATH);

  if (!Array.isArray(colorings)) {
    return document;
  }

  return setAt(
    document,
    COLORINGS_PATH,
    colorings.map((coloring) => {
      const key = isJsonObject(coloring) ? coloring["key"] : undefined;
      const scale = isString(key) ? scales.get(key) : undefined;

      return isJsonObject(coloring) && scale !== undefined
        ? setAt(
            coloring,
            SCALE_PATH,
            scale.map(([state, color]) => [state, color]),
          )
        : coloring;
    }),
  );
}
