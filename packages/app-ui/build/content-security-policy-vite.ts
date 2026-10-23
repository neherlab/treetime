import type { Plugin } from "vite";

export interface ContentSecurityPolicyOptions {
  webAssembly?: boolean;
}

export type ContentSecurityPolicyDirective = readonly [name: string, sources: readonly string[]];

export function contentSecurityPolicyDirectives(
  devServer: boolean,
  { webAssembly = false }: ContentSecurityPolicyOptions = {},
): ContentSecurityPolicyDirective[] {
  return [
    ["default-src", ["'self'"]],
    [
      "script-src",
      ["'self'", ...(devServer ? ["'unsafe-inline'"] : []), ...(webAssembly ? ["'wasm-unsafe-eval'"] : [])],
    ],
    ["style-src", ["'self'", "'unsafe-inline'"]],
    ["img-src", ["'self'", "data:", "blob:"]],
    ["font-src", ["'self'", "data:"]],
    ["connect-src", ["'self'"]],
    ["worker-src", ["'self'", "blob:"]],
    ["object-src", ["'none'"]],
    ["base-uri", ["'none'"]],
    ["form-action", ["'none'"]],
  ];
}

export function contentSecurityPolicy(devServer: boolean, options: ContentSecurityPolicyOptions = {}): string {
  return contentSecurityPolicyDirectives(devServer, options)
    .map(([name, sources]) => `${name} ${sources.join(" ")}`)
    .join("; ");
}

export function contentSecurityPolicyMeta(options: ContentSecurityPolicyOptions = {}): Plugin {
  let devServer = false;

  return {
    name: "treetime-content-security-policy",
    configResolved(config) {
      devServer = config.command === "serve";
    },
    transformIndexHtml() {
      return [
        {
          tag: "meta",
          attrs: { "http-equiv": "Content-Security-Policy", content: contentSecurityPolicy(devServer, options) },
          injectTo: "head-prepend",
        },
      ];
    },
  };
}
