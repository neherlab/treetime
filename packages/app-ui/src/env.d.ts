/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly VITE_TREETIME_DOCS_URL?: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}

declare module "*.css" {}
