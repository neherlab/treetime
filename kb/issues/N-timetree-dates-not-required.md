# --dates not required, misleading error when omitted

`--dates` is optional in clap, but the timetree command requires date constraints to function. When omitted, the error is "No variation in sampling dates! Please specify your clock rate explicitly." which does not mention that `--dates` was not provided.

The flag should be required for the timetree command, or the error message should suggest providing `--dates`.
