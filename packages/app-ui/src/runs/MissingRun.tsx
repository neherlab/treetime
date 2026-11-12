import { Button } from "@neherlab/design/ui/button";
import { Empty, EmptyContent, EmptyDescription, EmptyHeader, EmptyTitle } from "@neherlab/design/ui/empty";
import { Link } from "@tanstack/react-router";
import { useCallback } from "react";
import { useDropzone } from "react-dropzone";

import { dropzoneAccept } from "../analysis/fileAccept";
import { useHost } from "../host-context";
import { useImportRun } from "../start/useImportRun";

const PACKAGE_EXTENSIONS = ["zip"];

export function MissingRun() {
  const desktop = useHost() !== null;
  const importRun = useImportRun();

  const onDropAccepted = useCallback(
    (files: File[]) => {
      const [file] = files;

      if (file !== undefined) {
        void importRun(file);
      }
    },
    [importRun],
  );

  const { getInputProps, open: choosePackage } = useDropzone({
    noClick: true,
    noDrag: true,
    noKeyboard: true,
    multiple: false,
    onDropAccepted,
    ...dropzoneAccept(PACKAGE_EXTENSIONS),
  });

  return (
    <Empty className="py-10">
      <EmptyHeader>
        <EmptyTitle>{desktop ? "This run is not in the runs folder" : "This run is not in this tab"}</EmptyTitle>
        <EmptyDescription>
          {desktop
            ? "The run was deleted, or the runs folder changed since the run."
            : "The tab was reloaded or closed since the run, or the run was deleted."}
        </EmptyDescription>
      </EmptyHeader>
      <EmptyContent className="flex-row justify-center">
        <input {...getInputProps()} aria-label="Run package to import" />
        <Button variant="outline" nativeButton={false} render={<Link to="/" />}>
          Open the start page
        </Button>
        <Button type="button" variant="outline" onClick={choosePackage}>
          Import a run package
        </Button>
      </EmptyContent>
    </Empty>
  );
}
