import { errorMessage, type IndexedDataset } from "@neherlab/app-contracts";
import { Alert, AlertDescription } from "@neherlab/design/ui/alert";
import { Button } from "@neherlab/design/ui/button";
import { Spinner } from "@neherlab/design/ui/spinner";
import { useToastManager } from "@neherlab/design/ui/toast";
import { useCallback, useMemo, useState } from "react";

import type { ExampleIndex } from "../analysis/exampleData";
import { ExamplesDownload } from "../analysis/ExamplesDownload";
import { useExampleIndex, useOpenExample } from "../analysis/useExampleData";
import { yearDateText } from "../format";
import { galleryOrganisms } from "./gallery";

export function ExampleGallery() {
  const { data: examples, error } = useExampleIndex();
  const organisms = useMemo(() => galleryOrganisms(examples?.index?.datasets ?? []), [examples]);

  return (
    <section aria-labelledby="example-gallery-title" className="grid gap-4">
      <div className="grid gap-1">
        <h2 id="example-gallery-title" className="text-lg font-bold">
          Example datasets
        </h2>
        <p className="text-muted-foreground text-sm">
          Each dataset opens as a new time tree analysis with its tree, alignment and metadata.
        </p>
      </div>
      {error !== null && (
        <Alert variant="destructive">
          <AlertDescription>The examples cannot be listed: {error.message}</AlertDescription>
        </Alert>
      )}
      <ExamplesDownload />
      {examples !== undefined &&
        organisms.map((organism) => (
          <div key={organism.label} className="grid gap-2">
            <h3 className="text-sm font-bold">{organism.label}</h3>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(13rem,1fr))] gap-2">
              {organism.datasets.map((dataset) => (
                <DatasetCard key={dataset.name} examples={examples} dataset={dataset} />
              ))}
            </ul>
          </div>
        ))}
    </section>
  );
}

function DatasetCard({ examples, dataset }: { examples: ExampleIndex; dataset: IndexedDataset }) {
  const openExample = useOpenExample();
  const toasts = useToastManager();
  const [opening, setOpening] = useState(false);

  const open = useCallback(async () => {
    setOpening(true);

    try {
      await openExample(examples, dataset);
    } catch (error: unknown) {
      toasts.add({ title: `${dataset.name} cannot be opened`, description: errorMessage(error) });
    } finally {
      setOpening(false);
    }
  }, [dataset, examples, openExample, toasts]);

  const onOpen = useCallback(() => void open(), [open]);

  return (
    <li className="bg-card grid gap-2 rounded-lg border p-3" aria-label={dataset.name}>
      <div className="flex items-baseline justify-between gap-2">
        <span className="font-mono text-xs">{dataset.name}</span>
        <Button type="button" variant="outline" size="xs" disabled={opening} onClick={onOpen}>
          {opening && <Spinner />}
          Open
        </Button>
      </div>
      <dl className="text-muted-foreground grid grid-cols-[auto_1fr] gap-x-3 gap-y-0.5 text-xs">
        <dt>Samples</dt>
        <dd className="text-foreground tabular-nums">{dataset.samples.toLocaleString("en-US")}</dd>
        {dataset.dates !== undefined && (
          <>
            <dt>Dates</dt>
            <dd className="text-foreground tabular-nums">
              {yearDateText(dataset.dates.earliest)} to {yearDateText(dataset.dates.latest)}
            </dd>
          </>
        )}
        {dataset.sequence_length !== undefined && (
          <>
            <dt>Sequence length</dt>
            <dd className="text-foreground tabular-nums">{dataset.sequence_length.toLocaleString("en-US")} sites</dd>
          </>
        )}
      </dl>
    </li>
  );
}
