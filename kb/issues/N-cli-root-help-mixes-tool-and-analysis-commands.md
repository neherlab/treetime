# Root help mixes tool commands with analysis commands

The subcommand list of `treetime --help` shows the tool commands (`examples`, `schema`, `completions`, `help-markdown`) at the same level as the analysis commands (`timetree`, `ancestral`, `clock`, `mugration`, `homoplasy`, `optimize`, `prune`, `pipeline`), without a group heading or a label. A reader must know which commands analyze data and which serve setup and tooling. The order of the list follows the declaration order of `enum TreetimeCommands` in `packages/app-cli/src/cli/treetime_cli.rs`.

clap prints all subcommands under one `Commands:` heading and has no subcommand groups, so a grouped list needs a custom help template (`help_template` with the groups written out) or a hand-written subcommand section, and the Markdown reference of `treetime help-markdown` needs the same groups.

## Fix direction

Show two groups in the root help, analysis commands and tool commands, and test that every subcommand belongs to exactly one group.
