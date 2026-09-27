import { ThemeProvider as NextThemesProvider } from "next-themes";

const DATA_BLOCK_SCRIPT = { type: "application/json" } as const;

export function ThemeProvider({ children }: { children: React.ReactNode }) {
  return (
    <NextThemesProvider
      attribute="class"
      defaultTheme="system"
      disableTransitionOnChange
      scriptProps={DATA_BLOCK_SCRIPT}
    >
      {children}
    </NextThemesProvider>
  );
}
