import type { ReactNode } from "react";

import { cn } from "../ui/cn";

export function PageShell({ children, className }: { children: ReactNode; className?: string }) {
  return (
    <div
      className={cn("mx-auto grid w-full max-w-[110rem] grid-cols-1 content-start gap-4 px-5 pt-5 pb-16", className)}
    >
      {children}
    </div>
  );
}

export function PageHeading({
  title,
  description,
  actions,
}: {
  title: ReactNode;
  description?: ReactNode;
  actions?: ReactNode;
}) {
  return (
    <header className="flex flex-wrap items-start gap-4">
      <div className="grid min-w-0 gap-1">
        <h1 className="font-heading text-2xl leading-tight font-bold">{title}</h1>
        {description !== undefined && <p className="text-muted-foreground text-sm">{description}</p>}
      </div>
      {actions !== undefined && <div className="ml-auto flex flex-wrap items-center gap-1.5">{actions}</div>}
    </header>
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
