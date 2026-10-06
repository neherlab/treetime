import { decodeColorByGenotype, encodeColorByGenotype, isColorByGenotype } from "auspice/src/util/getGenotype";
import { nucleotide_gene } from "auspice/src/util/globals";

export function nucleotideColorBy(position: number): string | undefined {
  return encodeColorByGenotype({ gene: nucleotide_gene, positions: [position] }) ?? undefined;
}

export function nucleotidePosition(colorBy: string): number | undefined {
  if (!isColorByGenotype(colorBy)) {
    return undefined;
  }

  const genotype = decodeColorByGenotype(colorBy);

  if (genotype === null || genotype.gene !== nucleotide_gene || genotype.positions.length !== 1) {
    return undefined;
  }

  return genotype.positions[0];
}
