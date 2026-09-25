import * as z from "zod";

export interface TreeNode {
  name: string;
  children: TreeNode[];
  div: number | undefined;
  date: number | undefined;
  dateInterval: readonly [number, number] | undefined;
  excluded: boolean | undefined;
  traits: ReadonlyMap<string, TraitValue>;
  mutations: readonly string[];
}

interface TraitValue {
  value: string;
  confidence: Readonly<Record<string, number>> | undefined;
}

export interface ResultTree {
  root: TreeNode;
  nodes: readonly TreeNode[];
  tips: readonly TreeNode[];
  parents: ReadonlyMap<TreeNode, TreeNode>;
  colorings: readonly Coloring[];
  defaultColorBy: string | undefined;
}

interface Coloring {
  key: string;
  title: string;
  type: string;
}

export type AuspiceDocument = z.infer<typeof zAuspiceDocument>;

const RESERVED_ATTRS = new Set(["div", "num_date", "bad_branch"]);

export function parseAuspiceJson(text: string): AuspiceDocument {
  return zAuspiceDocument.parse(JSON.parse(text));
}

export function readAuspiceTree(document: AuspiceDocument): ResultTree {
  const root = toNode(document.tree);
  const nodes = preorder(root);
  const parents = new Map(nodes.flatMap((node) => node.children.map((child) => [child, node] as const)));

  return {
    root,
    nodes,
    tips: nodes.filter((node) => node.children.length === 0),
    parents,
    colorings: document.meta.colorings ?? [],
    defaultColorBy: document.meta.display_defaults?.color_by ?? undefined,
  };
}

export function tipNames(node: TreeNode): string[] {
  return node.children.length === 0 ? [node.name] : node.children.flatMap(tipNames);
}

export function initialColorBy(tree: ResultTree, preferred: readonly string[]): string | undefined {
  const keys = new Set(tree.colorings.map((coloring) => coloring.key));

  return preferred.find((key) => keys.has(key));
}

const zAttribute = z.looseObject({
  value: z.union([z.string(), z.number(), z.boolean()]),
  confidence: z.union([z.record(z.string(), z.number()), z.tuple([z.number(), z.number()])]).optional(),
});

const zRawNode = z.looseObject({
  name: z.string(),
  node_attrs: z.record(z.string(), z.union([z.number(), z.string(), zAttribute])).optional(),
  branch_attrs: z.looseObject({ mutations: z.record(z.string(), z.array(z.string())).optional() }).nullish(),
  get children() {
    return z.array(zRawNode).optional();
  },
});

type RawNode = z.infer<typeof zRawNode>;

const zNodeAttrs = z.looseObject({
  div: z.number().optional(),
  num_date: z.object({ value: z.number(), confidence: z.tuple([z.number(), z.number()]).optional() }).optional(),
  bad_branch: z.object({ value: z.string() }).optional(),
});

const zAuspiceDocument = z.looseObject({
  meta: z.looseObject({
    colorings: z.array(z.looseObject({ key: z.string(), title: z.string(), type: z.string() })).nullish(),
    display_defaults: z.looseObject({ color_by: z.string().nullish() }).nullish(),
  }),
  tree: zRawNode,
});

function toNode(raw: RawNode): TreeNode {
  const attrs = zNodeAttrs.parse(raw.node_attrs ?? {});

  return {
    name: raw.name,
    children: (raw.children ?? []).map(toNode),
    div: attrs.div,
    date: attrs.num_date?.value,
    dateInterval: attrs.num_date?.confidence,
    excluded: attrs.bad_branch === undefined ? undefined : attrs.bad_branch.value === "Yes",
    traits: new Map(
      Object.entries(raw.node_attrs ?? {})
        .filter(([key]) => !RESERVED_ATTRS.has(key))
        .flatMap(([key, value]) => {
          const trait = zAttribute.safeParse(value);

          return trait.success ? [[key, traitValue(trait.data)] as const] : [];
        }),
    ),
    mutations: raw.branch_attrs?.mutations?.["nuc"] ?? [],
  };
}

function traitValue(attribute: z.infer<typeof zAttribute>): TraitValue {
  const confidence = attribute.confidence;

  return { value: String(attribute.value), confidence: Array.isArray(confidence) ? undefined : confidence };
}

function preorder(root: TreeNode): TreeNode[] {
  return [root, ...root.children.flatMap(preorder)];
}
