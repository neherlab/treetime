import type { ImportDeclaration, Program } from "estree";
import type { Root } from "mdast";
import type { MdxJsxAttribute, MdxJsxFlowElement, MdxJsxTextElement } from "mdast-util-mdx-jsx";
import type { MdxjsEsm } from "mdast-util-mdxjs-esm";
import { visit } from "unist-util-visit";

import { type ProvidedImport, providedImport } from "./provided-imports";

export function remarkProvidedComponents() {
  return (tree: Root): void => {
    const imported = importedNames(tree);
    const used = new Map<string, ProvidedImport>();

    visit(tree, (node) => {
      if ((node.type !== "mdxJsxFlowElement" && node.type !== "mdxJsxTextElement") || node.name === null) {
        return;
      }

      const provided = imported.has(node.name) ? undefined : providedImport(node.name);

      if (provided !== undefined) {
        used.set(node.name, provided);
        addHydration(node, provided);
      }
    });

    tree.children.unshift(...[...used].map(([name, provided]) => importNode(name, provided)));
  };
}

function importedNames(tree: Root): Set<string> {
  return new Set(
    tree.children
      .flatMap((node) => (node.type === "mdxjsEsm" ? (node.data?.estree?.body ?? []) : []))
      .flatMap((statement) => (statement.type === "ImportDeclaration" ? statement.specifiers : []))
      .map((specifier) => specifier.local.name),
  );
}

function addHydration(node: MdxJsxFlowElement | MdxJsxTextElement, provided: ProvidedImport): void {
  const hydrated = node.attributes.some(
    (attribute) => attribute.type === "mdxJsxAttribute" && attribute.name.startsWith("client:"),
  );

  if (provided.client !== undefined && !hydrated) {
    const directive: MdxJsxAttribute = { type: "mdxJsxAttribute", name: `client:${provided.client}`, value: null };
    node.attributes.push(directive);
  }
}

function importNode(name: string, { source, exported }: ProvidedImport): MdxjsEsm {
  const local = { type: "Identifier", name } as const;

  const declaration: ImportDeclaration = {
    type: "ImportDeclaration",
    specifiers: [
      exported === "default"
        ? { type: "ImportDefaultSpecifier", local }
        : { type: "ImportSpecifier", imported: { type: "Identifier", name: exported }, local },
    ],
    source: { type: "Literal", value: source, raw: JSON.stringify(source) },
    attributes: [],
  };

  const program: Program = { type: "Program", sourceType: "module", body: [declaration] };
  const specifier = exported === "default" ? name : `{ ${exported} }`;

  return { type: "mdxjsEsm", value: `import ${specifier} from ${JSON.stringify(source)};`, data: { estree: program } };
}
