import type { IndexedDataset } from "@neherlab/app-contracts";

import { galleryOrganisms } from "../gallery";

function dataset(name: string, organism: string, samples: number, segment?: string): IndexedDataset {
  const indexed: IndexedDataset = { name, organism, samples, inputs: [{ kind: "tree", path: `${name}/tree.nwk` }] };

  if (segment !== undefined) {
    indexed.segment = segment;
  }

  return indexed;
}

describe("example gallery", () => {
  test("datasets are grouped by organism and segment, in the order organism, then size", () => {
    const datasets = [
      dataset("zika/86", "Zika virus", 86),
      dataset("ebola/362", "Ebola virus", 362),
      dataset("flu/h3n2/ha/20", "Influenza A virus H3N2", 19, "HA"),
      dataset("ebola/20", "Ebola virus", 20),
      dataset("flu/h3n2/na/20", "Influenza A virus H3N2", 20, "NA"),
      dataset("zika/20", "Zika virus", 20),
      dataset("ebola/100", "Ebola virus", 100),
    ];

    expect(
      galleryOrganisms(datasets).map((group) => [group.label, group.datasets.map((member) => member.name)]),
    ).toStrictEqual([
      ["Ebola virus", ["ebola/20", "ebola/100", "ebola/362"]],
      ["Influenza A virus H3N2, segment HA", ["flu/h3n2/ha/20"]],
      ["Influenza A virus H3N2, segment NA", ["flu/h3n2/na/20"]],
      ["Zika virus", ["zika/20", "zika/86"]],
    ]);
  });
});
