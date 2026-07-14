# Validate finite and representable grid spacing

Make every public grid-producing boundary enforce the documented numeric-domain invariants before constructing a value.

## Acceptance criteria

- Reject NaN and positive or negative infinity for origins, endpoints, and spacing.
- Reject zero and negative spacing.
- Reject ranges whose generated adjacent coordinates are not strictly increasing in the represented type.
- Replace fallible numeric-conversion `unwrap()` calls with contextual errors.
- Validate derived deserialization or replace it with a validated implementation so serialized input cannot bypass the invariant.
- Validate `GridFn::from_arrays_nonuniform()` before any float-to-`usize` conversion.
- Cover every constructor, deserialization, and `GridFn` boundary with NaN, infinity, reversed range, zero spacing, extreme finite span, and repeated-representable-coordinate cases.

## Related issues

Source: [kb/issues/M-grid-invalid-numeric-domain.md](../issues/M-grid-invalid-numeric-domain.md)
