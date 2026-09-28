import type { DropzoneOptions } from "react-dropzone";

const ANY_APPLICATION_TYPE = "application/*";

export function dropzoneAccept(extensions: readonly string[]): Pick<DropzoneOptions, "accept"> {
  return extensions.length === 0
    ? {}
    : { accept: { [ANY_APPLICATION_TYPE]: extensions.map((extension) => `.${extension}`) } };
}
