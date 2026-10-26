import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@neherlab/app-ui/ui/table";
import type { ComponentProps } from "react";

import { CalloutContainer, CalloutDescription, CalloutTitle } from "./Callout";

export const mdxComponents = {
  CalloutContainer,
  CalloutTitle,
  CalloutDescription,
  table: DocsTable,
  thead: TableHeader,
  tbody: TableBody,
  tr: TableRow,
  th: TableHead,
  td: TableCell,
};

function DocsTable(props: ComponentProps<"table">) {
  return (
    <div className="not-prose my-6">
      <Table {...props} />
    </div>
  );
}
