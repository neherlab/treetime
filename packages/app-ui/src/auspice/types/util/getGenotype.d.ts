export declare function isColorByGenotype(colorBy: string): boolean;

export declare function encodeColorByGenotype(genotype: { gene?: string; positions: readonly number[] }): string | null;

export declare function decodeColorByGenotype(colorBy: string): DecodedGenotype | null;

interface DecodedGenotype {
  gene: string;
  positions: number[];
  aa: boolean;
}
