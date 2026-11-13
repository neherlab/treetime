import type { CheckConfigResponse, SettingDifference } from "@neherlab/app-contracts";

export function validCheck(differences: SettingDifference[]): CheckConfigResponse {
  return {
    status: "valid",
    command: "timetree",
    config: {},
    code: {
      command_line: [],
      command_line_text: "",
      yaml: [],
      yaml_text: "",
      config_file: "",
      config_command: "",
      methods: "",
      citation: { text: "", short: "", doi: "", url: "" },
    },
    checks: [],
    choices: [],
    effective: [],
    command_lines: [],
    run_differences: differences,
  };
}

export function invalidCheck(): CheckConfigResponse {
  return {
    status: "invalid",
    command: "timetree",
    message: "bad",
    causes: [],
    problems: [],
    messages: [],
    checks: [],
    choices: [],
    effective: [],
    command_lines: [],
  };
}
