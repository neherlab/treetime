import { tipNames, type ResultTree, type TreeNode } from "./tree";

export interface CladeIndex {
  keyOf: ReadonlyMap<TreeNode, string>;
  nodeOf: ReadonlyMap<string, TreeNode>;
}

export interface MatchedClade {
  key: string;
  first: TreeNode;
  second: TreeNode;
}

function cladeKey(names: readonly string[]): string {
  return JSON.stringify(names.toSorted());
}

export function indexClades(tree: ResultTree): CladeIndex {
  const keyOf = new Map(tree.nodes.map((node) => [node, cladeKey(tipNames(node))] as const));
  const nodeOf = new Map<string, TreeNode>();

  for (const [node, key] of keyOf) {
    if (!nodeOf.has(key)) {
      nodeOf.set(key, node);
    }
  }

  return { keyOf, nodeOf };
}

export function matchAncestors(first: CladeIndex, second: CladeIndex): MatchedClade[] {
  return [...first.nodeOf].flatMap(([key, node]) => {
    const other = second.nodeOf.get(key);

    return node.children.length > 0 && other !== undefined && other.children.length > 0
      ? [{ key, first: node, second: other }]
      : [];
  });
}
