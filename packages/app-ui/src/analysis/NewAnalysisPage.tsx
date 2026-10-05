import { useDraftStore } from "../store/draft";
import { DraftForm } from "./DraftForm";

export function NewAnalysisPage() {
  const command = useDraftStore((state) => state.draft.command);
  const epoch = useDraftStore((state) => state.epoch);

  return <DraftForm key={`${command}:${epoch}`} command={command} />;
}
