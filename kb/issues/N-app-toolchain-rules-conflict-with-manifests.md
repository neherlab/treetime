# Application toolchain rules conflict with package manifests

Project rules require TypeScript 5.9 and React 18, while application manifests pin TypeScript 6.0.3 and React 19.2.6. Both sources are explicit; changing either affects the governing development contract.

## Options

- O1. Align package manifests and lockfiles to TypeScript 5.9 and React 18.
- O2. Amend the governing rules through the project decision process to authorize the manifest versions.

No ticket is ready until the project chooses which source is authoritative.

## Recommendation

Align manifests to the currently governing TypeScript 5.9 and React 18 rules unless the project explicitly approves a rule change.

## Related issues

- [H-core-multi-client-architecture-library-purity.md](H-core-multi-client-architecture-library-purity.md)
