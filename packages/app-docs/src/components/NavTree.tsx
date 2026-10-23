import { LocaleProvider, m, type Locale } from "@neherlab/app-ui/i18n";
import {
  SidebarGroup,
  SidebarGroupLabel,
  SidebarMenu,
  SidebarMenuButton,
  SidebarMenuItem,
} from "@neherlab/app-ui/ui/sidebar";
import { match } from "ts-pattern";
import ChevronRight from "~icons/lucide/chevron-right";

import { containsUrl, type NavNode } from "../lib/nav";

export function NavTree({ nodes, currentUrl, locale }: { nodes: NavNode[]; currentUrl: string; locale: Locale }) {
  return (
    <LocaleProvider locale={locale}>
      <nav aria-label={m.docs_navigation({}, { locale })}>
        <NavList nodes={nodes} currentUrl={currentUrl} />
      </nav>
    </LocaleProvider>
  );
}

function NavList({ nodes, currentUrl }: { nodes: NavNode[]; currentUrl: string }) {
  return (
    <SidebarMenu>
      {nodes.map((node) => (
        <SidebarMenuItem key={`${node.kind}:${node.name}`}>
          <NavEntry node={node} currentUrl={currentUrl} />
        </SidebarMenuItem>
      ))}
    </SidebarMenu>
  );
}

function NavEntry({ node, currentUrl }: { node: NavNode; currentUrl: string }) {
  return match(node)
    .with({ kind: "page" }, (page) => <PageLink name={page.name} url={page.url} currentUrl={currentUrl} />)
    .with({ kind: "separator" }, (separator) => (
      <p className="text-muted-foreground mt-3 px-2 text-xs font-bold">{separator.name}</p>
    ))
    .with({ kind: "folder" }, (folder) => <NavFolder folder={folder} currentUrl={currentUrl} />)
    .exhaustive();
}

function NavFolder({ folder, currentUrl }: { folder: Extract<NavNode, { kind: "folder" }>; currentUrl: string }) {
  return (
    <SidebarGroup className="p-0">
      <details open={containsUrl(folder, currentUrl)}>
        <SidebarGroupLabel
          // oxlint-disable-next-line jsx-a11y/control-has-associated-label -- the summary gets its text from the label children
          render={<summary />}
          className="hover:bg-sidebar-accent cursor-pointer list-none gap-1 [&::-webkit-details-marker]:hidden"
        >
          <ChevronRight aria-hidden className="transition-transform [details[open]>summary>&]:rotate-90" />
          {folder.name}
        </SidebarGroupLabel>
        <div className="border-sidebar-border ml-3.5 border-l pl-2">
          {folder.url !== undefined && <PageLink name={folder.name} url={folder.url} currentUrl={currentUrl} />}
          <NavList nodes={folder.children} currentUrl={currentUrl} />
        </div>
      </details>
    </SidebarGroup>
  );
}

function PageLink({ name, url, currentUrl }: { name: string; url: string; currentUrl: string }) {
  const current = url === currentUrl;

  return (
    <SidebarMenuButton
      isActive={current}
      // oxlint-disable-next-line jsx-a11y/control-has-associated-label -- the link gets its text from the button children
      render={<a href={url} aria-current={current ? "page" : undefined} />}
    >
      <span>{name}</span>
    </SidebarMenuButton>
  );
}
