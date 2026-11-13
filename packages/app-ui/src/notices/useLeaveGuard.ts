import { useEffect } from "react";

import { useHost } from "../host-context";
import { useListedRuns } from "../shell/useListedRuns";
import { useDraftStore } from "../store/draft";
import { formUploads, leaveLosesWork, useTabWorkStore } from "./leaveGuard";

export function useLeaveGuard(): void {
  const host = useHost();
  const { all } = useListedRuns();
  const imports = useTabWorkStore((state) => state.imports);
  const setHoldsWork = useTabWorkStore((state) => state.setHoldsWork);
  const uploads = useDraftStore((state) => formUploads(state.draft));

  const holdsWork =
    host === null &&
    leaveLosesWork({
      runs: all.length,
      imports,
      formUploads: uploads,
    });

  useEffect(() => {
    setHoldsWork(holdsWork);
  }, [holdsWork, setHoldsWork]);

  useEffect(() => {
    if (!holdsWork) {
      return undefined;
    }

    const warn = (event: BeforeUnloadEvent) => {
      event.preventDefault();
    };

    window.addEventListener("beforeunload", warn);

    return () => {
      window.removeEventListener("beforeunload", warn);
    };
  }, [holdsWork]);
}
