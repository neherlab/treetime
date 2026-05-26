/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly TREETIME_DEBUG_FETCH?: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}

declare module "*.css" {}
