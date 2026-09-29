import type { Plugin } from "vite";

export function contentSecurityPolicy(devServer: boolean): string {
  const scriptSources = devServer ? "'self' 'unsafe-inline'" : "'self'";

  return [
    "default-src 'self'",
    `script-src ${scriptSources}`,
    "style-src 'self' 'unsafe-inline'",
    "img-src 'self' data: blob:",
    "font-src 'self' data:",
    "connect-src 'self'",
    "worker-src 'self' blob:",
    "object-src 'none'",
    "base-uri 'none'",
    "form-action 'none'",
  ].join("; ");
}

export function contentSecurityPolicyMeta(): Plugin {
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
          attrs: { "http-equiv": "Content-Security-Policy", content: contentSecurityPolicy(devServer) },
          injectTo: "head-prepend",
        },
      ];
    },
  };
}
