import { Meter as MeterPrimitive } from "@base-ui/react/meter";
import { cn } from "cn";

function Meter({ className, children, ...props }: MeterPrimitive.Root.Props) {
  return (
    <MeterPrimitive.Root data-slot="meter" className={cn("flex flex-wrap gap-x-3 gap-y-1.5", className)} {...props}>
      {children}
      <MeterTrack>
        <MeterIndicator />
      </MeterTrack>
    </MeterPrimitive.Root>
  );
}

function MeterTrack({ className, ...props }: MeterPrimitive.Track.Props) {
  return (
    <MeterPrimitive.Track
      data-slot="meter-track"
      className={cn("bg-muted relative h-1.5 w-full overflow-hidden rounded-full", className)}
      {...props}
    />
  );
}

function MeterIndicator({ className, ...props }: MeterPrimitive.Indicator.Props) {
  return (
    <MeterPrimitive.Indicator
      data-slot="meter-indicator"
      className={cn("bg-muted-foreground h-full transition-all", className)}
      {...props}
    />
  );
}

function MeterLabel({ className, ...props }: MeterPrimitive.Label.Props) {
  return <MeterPrimitive.Label data-slot="meter-label" className={cn("text-sm font-medium", className)} {...props} />;
}

function MeterValue({ className, ...props }: MeterPrimitive.Value.Props) {
  return (
    <MeterPrimitive.Value
      data-slot="meter-value"
      className={cn("text-muted-foreground ml-auto text-sm tabular-nums", className)}
      {...props}
    />
  );
}

export { Meter, MeterLabel, MeterValue };
