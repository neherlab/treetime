import { fileURLToPath } from "node:url";

import Icons from "unplugin-icons/vite";
import type { PluginOption } from "vite";

const ICON_SIZE = "24";

export function icons(): PluginOption {
  return Icons({
    compiler: "jsx",
    jsx: "react",
    collectionsNodeResolvePath: fileURLToPath(new URL("..", import.meta.url)),
    iconCustomizer(_collection, _icon, props) {
      props["width"] = ICON_SIZE;
      props["height"] = ICON_SIZE;
    },
  });
}
