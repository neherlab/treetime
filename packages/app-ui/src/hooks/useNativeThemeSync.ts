import { zUiTheme } from "@neherlab/app-contracts";
import { useTheme } from "next-themes";
import { useEffect } from "react";

import { useHost } from "../host-context";

export function useNativeThemeSync() {
  const host = useHost();
  const { theme } = useTheme();
  const native = zUiTheme.options.map((option) => option.value).find((value) => value === theme);

  useEffect(() => {
    if (host !== null && native !== undefined) {
      void host.setNativeTheme(native);
    }
  }, [host, native]);
}
