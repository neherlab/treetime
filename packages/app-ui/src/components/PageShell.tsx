import { cn } from "cn";
import type { ReactNode } from "react";

export function PageShell({ children, className }: { children: ReactNode; className?: string }) {
  return (
    <div
      className={cn("mx-auto grid w-full max-w-[110rem] grid-cols-1 content-start gap-4 px-5 pt-5 pb-16", className)}
    >
      {children}
    </div>
  );
}

export function LoadingState({ text }: { text: string }) {
  return (
    <output className="text-muted-foreground flex flex-col items-center justify-center gap-3 p-10 text-sm">
      <img src={`${import.meta.env.BASE_URL}logo-loading.svg`} alt="" className="size-12" />
      {text}
    </output>
  );
}
