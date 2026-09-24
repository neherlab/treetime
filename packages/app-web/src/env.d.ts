/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly VITE_TREETIME_DEBUG_FETCH?: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}

declare module "*.css" {}
