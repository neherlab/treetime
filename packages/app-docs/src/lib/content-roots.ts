export const CONTENT_ROOTS = {
  docs: "../../docs/",
  generated: "src/generated/",
} as const;

export type ContentRoot = keyof typeof CONTENT_ROOTS;

export const EDIT_ROOT = "https://github.com/neherlab/treetime/edit/rust/docs/";
