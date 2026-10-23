import { readdir, readFile, stat, writeFile } from "node:fs/promises";
import { join, relative } from "node:path";
import { constants, gzipSync } from "node:zlib";

const COMPRESSIBLE = /\.(?:html|css|js|json|svg|xml)$/u;

const THRESHOLD = 1024;

export async function writeCompressedCopies(root: string, skippedFolders: readonly string[]): Promise<number> {
  const files = await compressibleFiles(root, skippedFolders);

  await Promise.all(
    files.map(async (file) => {
      const content = await readFile(file);

      await writeFile(`${file}.gz`, gzipSync(content, { level: constants.Z_BEST_COMPRESSION }));
    }),
  );

  return files.length;
}

async function compressibleFiles(root: string, skippedFolders: readonly string[]): Promise<string[]> {
  const entries = await readdir(root, { recursive: true, withFileTypes: true });

  const candidates = entries.flatMap((entry) => {
    const file = join(entry.parentPath, entry.name);
    const skipped = skippedFolders.some((folder) => relative(root, file).startsWith(`${folder}/`));

    return entry.isFile() && COMPRESSIBLE.test(entry.name) && !skipped ? [file] : [];
  });

  const sizes = await Promise.all(candidates.map(async (file) => (await stat(file)).size));

  return candidates.filter((_, index) => (sizes[index] ?? 0) > THRESHOLD);
}
