import type { ResultTree, TreeNode } from "./tree";

export const UNCERTAIN_STATE_PROBABILITY = 0.8;

export interface BranchMutations {
  name: string;
  tips: number;
  mutations: readonly string[];
}

export interface RecurrentSite {
  position: number;
  branches: number;
}

export interface StateChange {
  from: string;
  to: string;
  branches: number;
}

export interface AncestorState {
  name: string;
  tips: number;
  state: string;
  probability: number;
}

const MUTATION_POSITION = /^[^\d]*(?<position>\d+)[^\d]*$/u;

export function branchesByMutationCount(tree: ResultTree): BranchMutations[] {
  return tree.nodes
    .filter((node) => node.mutations.length > 0)
    .map((node) => ({ name: node.name, tips: tipCount(node), mutations: node.mutations }))
    .toSorted((left, right) => right.mutations.length - left.mutations.length || left.name.localeCompare(right.name));
}

export function recurrentSites(tree: ResultTree): RecurrentSite[] {
  const branchesPerSite = new Map<number, number>();

  for (const node of tree.nodes) {
    for (const position of new Set(node.mutations.map(mutationPosition))) {
      branchesPerSite.set(position, (branchesPerSite.get(position) ?? 0) + 1);
    }
  }

  return [...branchesPerSite]
    .flatMap(([position, branches]) => (branches > 1 ? [{ position, branches }] : []))
    .toSorted((left, right) => right.branches - left.branches || left.position - right.position);
}

function mutationPosition(mutation: string): number {
  const position = MUTATION_POSITION.exec(mutation)?.groups?.["position"];

  if (position === undefined) {
    throw new Error(`"${mutation}" does not name a sequence position`);
  }

  return Number(position);
}

export function stateChanges(tree: ResultTree, attribute: string): StateChange[] {
  const counts = new Map<string, StateChange>();

  for (const [child, parent] of tree.parents) {
    const from = parent.traits.get(attribute)?.value;
    const to = child.traits.get(attribute)?.value;

    if (from !== undefined && to !== undefined && from !== to) {
      const key = JSON.stringify([from, to]);
      const current = counts.get(key);
      counts.set(key, { from, to, branches: (current?.branches ?? 0) + 1 });
    }
  }

  return [...counts.values()].toSorted(
    (left, right) =>
      right.branches - left.branches || left.from.localeCompare(right.from) || left.to.localeCompare(right.to),
  );
}

export function ancestorStates(tree: ResultTree, attribute: string): AncestorState[] {
  return tree.nodes.flatMap((node) => {
    const trait = node.traits.get(attribute);
    const probability = trait?.confidence?.[trait.value];

    return node.children.length > 0 && trait !== undefined && probability !== undefined
      ? [{ name: node.name, tips: tipCount(node), state: trait.value, probability }]
      : [];
  });
}

export function uncertainAncestors(tree: ResultTree, attribute: string): AncestorState[] {
  return ancestorStates(tree, attribute)
    .filter((ancestor) => ancestor.probability < UNCERTAIN_STATE_PROBABILITY)
    .toSorted((left, right) => left.probability - right.probability);
}

export function tipCount(node: TreeNode): number {
  return node.children.length === 0 ? 1 : node.children.reduce((sum, child) => sum + tipCount(child), 0);
}
