import type { UploadAllowance } from "@neherlab/app-contracts";

import type { StoredFile } from "./fileSelection";

export interface SizedFile {
  name: string;
  size: number;
}

export async function storeWithinAllowance<T extends SizedFile>(
  files: readonly T[],
  allowance: () => Promise<UploadAllowance>,
  store: (files: readonly T[]) => Promise<StoredFile[]>,
): Promise<StoredFile[]> {
  const refusal = oversizedUpload(files, await allowance());

  if (refusal !== undefined) {
    throw new Error(refusal);
  }

  return store(files);
}

export function oversizedUpload(
  files: readonly SizedFile[],
  { remaining, limit_text }: UploadAllowance,
): string | undefined {
  let total = 0;

  for (const { name, size } of files) {
    total += size;

    if (total > remaining) {
      return `the inputs of a run are limited to ${limit_text} in total; \`${name}\` does not fit`;
    }
  }

  return undefined;
}
