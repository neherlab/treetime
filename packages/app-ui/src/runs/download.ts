export function saveBytes(bytes: Uint8Array, name: string, type: string): void {
  const url = URL.createObjectURL(new Blob([new Uint8Array(bytes)], { type }));
  const anchor = document.createElement("a");

  anchor.href = url;
  anchor.download = name;
  document.body.append(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
}
