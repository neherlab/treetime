import type { ReactNode } from "react";

import { Card } from "../ui/card";
import { ActionDemo, ChoiceDemo, FieldDemo, SelectDemo } from "./controls";
import { NavigationDemo, NoticeDemo, OverlayDemo, ProgressDemo, ReadoutDemo, StateDemo, TableDemo } from "./displays";

export function ComponentSection() {
  return (
    <div className="flex flex-col gap-4">
      <p className="text-muted-foreground max-w-prose text-sm">
        The application register. Each control is a component of the design package on Base UI or React Aria Components,
        so focus, keyboard use, and screen reader labels come from the library.
      </p>
      <div className="grid gap-4 @3xl:grid-cols-2 @6xl:grid-cols-3">
        <Demo title="Buttons">
          <ActionDemo />
        </Demo>
        <Demo title="Text fields">
          <FieldDemo />
        </Demo>
        <Demo title="Selects">
          <SelectDemo />
        </Demo>
        <Demo title="Choices">
          <ChoiceDemo />
        </Demo>
        <Demo title="Tabs, badges, states">
          <StateDemo />
        </Demo>
        <Demo title="Notices">
          <NoticeDemo />
        </Demo>
        <Demo title="Progress">
          <ProgressDemo />
        </Demo>
        <Demo title="Readouts">
          <ReadoutDemo />
        </Demo>
        <Demo title="Overlays">
          <OverlayDemo />
        </Demo>
        <Demo title="Navigation">
          <NavigationDemo />
        </Demo>
        <Demo title="Table" wide>
          <TableDemo />
        </Demo>
      </div>
    </div>
  );
}

function Demo({ title, wide = false, children }: { title: string; wide?: boolean; children: ReactNode }) {
  return (
    <Card size="sm" className={wide ? "gap-3 p-4 @3xl:col-span-2 @6xl:col-span-3" : "gap-3 p-4"}>
      <h3 className="text-muted-foreground text-xs font-semibold">{title}</h3>
      {children}
    </Card>
  );
}
