import { FileInputPanel } from "./FileInputPanel";
import { ParamForm } from "./ParamForm";
import { RunButton } from "./RunButton";

export function InputPanel() {
  return (
    <div className="flex flex-col gap-5 overflow-y-auto p-4">
      <FileInputPanel />
      <ParamForm />
      <RunButton />
    </div>
  );
}
