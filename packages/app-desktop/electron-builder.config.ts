import type { Configuration } from "electron-builder";

import desktop from "./package.json" with { type: "json" };

export default {
  appId: "org.neherlab.treetime",
  productName: "TreeTime",
  electronVersion: desktop.devDependencies.electron,
  artifactName: "treetime-desktop-${env.TREETIME_DESKTOP_TARGET}.${ext}",
  directories: {
    output: "../package",
  },
  asarUnpack: ["**/*.node"],
  npmRebuild: false,
  publish: null,
  electronFuses: {
    runAsNode: false,
    enableNodeOptionsEnvironmentVariable: false,
    enableNodeCliInspectArguments: false,
    enableCookieEncryption: true,
    enableEmbeddedAsarIntegrityValidation: true,
    onlyLoadAppFromAsar: true,
  },
  toolsets: {
    appimage: "1.0.3",
  },
  linux: {
    target: "AppImage",
    category: "Science",
    syncDesktopName: true,
  },
  win: {
    target: "zip",
  },
  mac: {
    target: "dmg",
    category: "public.app-category.education",
    identity: "-",
    hardenedRuntime: false,
  },
} satisfies Configuration;
