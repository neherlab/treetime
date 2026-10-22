# CLI flags parsed but ignored

No command has a known flag that it parses and never reads. Nothing prevents a new one: the command's raw argument struct (`Treetime<Command>ArgsRaw` in `packages/app-commands/src/commands/<command>/args.rs`) derives the deser `Serialize` and the clap parser, which read every field, so the compiler does not report a field that the run never reads. A user who sets such a flag gets no error and no effect, and the command reports success while it ignores part of its input.

## Guard against new unread flags

No test checks that a parsed argument field is used. Add one check that fails when a field of a command argument struct is neither read by the command nor rejected during validation, so that a new flag cannot be added without an effect.

## Related issues

- [M-cli-help-text-defects.md](M-cli-help-text-defects.md): help text defects of the command line
