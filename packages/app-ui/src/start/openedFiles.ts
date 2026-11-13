import type { DetectedFile, FileHead, FileRole } from "@neherlab/app-contracts";

const HEAD_BYTES = 65_536;

export interface RoledFile<F> {
  file: F;
  index: number;
  role: FileRole;
}

export async function fileHead(file: Blob): Promise<FileHead> {
  return { head: await file.slice(0, HEAD_BYTES).text(), complete: file.size <= HEAD_BYTES };
}

export function withRoles<F>(files: readonly F[], roles: readonly FileRole[]): Array<RoledFile<F>> {
  return files.map((file, index) => ({ file, index, role: roles[index] ?? "input" }));
}

export function detectedFiles(
  files: ReadonlyArray<RoledFile<{ name: string }>>,
  inputs: readonly DetectedFile[],
): DetectedFile[] {
  let next = 0;

  return files.map(({ file, role }): DetectedFile => {
    if (role !== "input") {
      return { path: file.name, kind: role, commands: [] };
    }

    const detected = inputs[next];
    next += 1;

    return detected ?? { path: file.name, kind: "unknown", commands: [], reason: "TreeTime did not read the file" };
  });
}
