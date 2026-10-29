import { cn } from "cn";
import { useCallback, useId, useState, type ReactNode } from "react";

import { Label } from "../ui/label";
import { Switch } from "../ui/switch";
import { ColorSection } from "./ColorSection";
import { ComponentSection } from "./ComponentSection";
import { IdentitySection } from "./IdentitySection";
import { RunMock } from "./RunMock";
import { SiteMock } from "./SiteMock";
import { TypeSection } from "./TypeSection";

export function DesignSpecimen() {
  const id = useId();
  const [both, setBoth] = useState(false);
  const toggle = useCallback((checked: boolean) => setBoth(checked), []);

  return (
    <div className="@container flex flex-col gap-12">
      <div className="bg-card sticky top-(--header-height) z-10 flex items-center gap-3 border-b py-3">
        <Switch id={id} checked={both} onCheckedChange={toggle} />
        <Label htmlFor={id}>Show the light and the dark theme side by side</Label>
      </div>
      <SpecimenSection title="Identity">
        <IdentitySection />
      </SpecimenSection>
      <SpecimenSection title="Color">
        <ThemeSamples both={both}>
          <ColorSection />
        </ThemeSamples>
      </SpecimenSection>
      <SpecimenSection title="Type">
        <ThemeSamples both={both}>
          <TypeSection />
        </ThemeSamples>
      </SpecimenSection>
      <SpecimenSection title="Components">
        <ThemeSamples both={both}>
          <ComponentSection />
        </ThemeSamples>
      </SpecimenSection>
      <SpecimenSection title="Application: a TreeTime run">
        <ThemeSamples both={both} stacked>
          <RunMock />
        </ThemeSamples>
      </SpecimenSection>
      <SpecimenSection title="Site: the lab home page">
        <ThemeSamples both={both} stacked>
          <SiteMock />
        </ThemeSamples>
      </SpecimenSection>
    </div>
  );
}

function SpecimenSection({ title, children }: { title: string; children: ReactNode }) {
  return (
    <section className="flex flex-col gap-5">
      <h2 className="font-serif text-3xl font-normal tracking-tight">{title}</h2>
      {children}
    </section>
  );
}

function ThemeSamples({ both, stacked = false, children }: { both: boolean; stacked?: boolean; children: ReactNode }) {
  if (!both) {
    return <div className="min-w-0">{children}</div>;
  }

  return (
    <div className={cn("grid gap-4", !stacked && "@5xl:grid-cols-2")}>
      <ThemeFrame scheme="light">{children}</ThemeFrame>
      <ThemeFrame scheme="dark">{children}</ThemeFrame>
    </div>
  );
}

function ThemeFrame({ scheme, children }: { scheme: "light" | "dark"; children: ReactNode }) {
  return (
    <section
      aria-label={scheme === "light" ? "Light theme" : "Dark theme"}
      className={cn(
        "bg-background text-foreground min-w-0 rounded-xl border p-4",
        scheme === "light" ? "scheme-light" : "scheme-dark",
      )}
    >
      <p className="text-muted-foreground mb-3 text-xs font-medium">{scheme === "light" ? "Light" : "Dark"}</p>
      {children}
    </section>
  );
}
