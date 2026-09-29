import { Link, useLocation } from "@tanstack/react-router";
import { useEffect } from "react";

import { Button } from "./ui/button";
import { Empty, EmptyContent, EmptyDescription, EmptyHeader, EmptyTitle } from "./ui/empty";

export function NotFoundPage() {
  const { pathname } = useLocation();

  useEffect(() => {
    console.error(`[TreeTime] no page matches the path "${pathname}"`);
  }, [pathname]);

  return (
    <Empty role="alert" className="py-10">
      <EmptyHeader>
        <EmptyTitle>Page not found</EmptyTitle>
        <EmptyDescription>No page matches the path {pathname}.</EmptyDescription>
      </EmptyHeader>
      <EmptyContent>
        <Button variant="outline" nativeButton={false} render={<Link to="/new" />}>
          New analysis
        </Button>
      </EmptyContent>
    </Empty>
  );
}
