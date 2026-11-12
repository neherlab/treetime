import type { IndexedDataset } from "@neherlab/app-contracts";

interface OrganismGroup {
  label: string;
  datasets: IndexedDataset[];
}

export function galleryOrganisms(datasets: readonly IndexedDataset[]): OrganismGroup[] {
  const sorted = datasets.toSorted(
    (left, right) => organismLabel(left).localeCompare(organismLabel(right)) || left.samples - right.samples,
  );

  return Object.entries(Object.groupBy(sorted, organismLabel)).flatMap(([label, members]) =>
    members === undefined ? [] : [{ label, datasets: members }],
  );
}

function organismLabel(dataset: IndexedDataset): string {
  return dataset.segment === undefined ? dataset.organism : `${dataset.organism}, segment ${dataset.segment}`;
}
