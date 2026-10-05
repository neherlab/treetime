import type { ActiveChoice, ChoiceOptionName, SettingChoice, SparseConfig } from "@neherlab/app-contracts";

import { sameJson } from "./json";

export function choiceRow(choice: SettingChoice, state: ChoiceState): ChoiceRow {
  const reported = state.reported?.find((active) => active.choice === choice.choice)?.option;
  const answered = state.checked !== undefined && sameJson(state.checked, state.config);
  const selected = (answered ? reported : state.picked) ?? reported ?? state.picked ?? choice.options[0]?.option;
  const option = choice.options.find((candidate) => candidate.option === selected);

  return {
    options: choice.options.map((candidate) => candidate.option),
    selected,
    shown: option?.settings ?? [],
  };
}

export interface ChoiceRow {
  options: ChoiceOptionName[];
  selected: ChoiceOptionName | undefined;
  shown: string[];
}

export interface ChoiceState {
  reported: readonly ActiveChoice[] | undefined;
  picked: ChoiceOptionName | undefined;
  checked: SparseConfig | undefined;
  config: SparseConfig;
}
