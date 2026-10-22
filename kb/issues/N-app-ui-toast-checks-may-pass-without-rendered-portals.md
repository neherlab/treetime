# Toast checks of the upload tests may pass without rendered portals

The upload tests check that no error toast appears with `queryByText(...)` equal to `null`. Toasts render through a Base UI portal. If happy-dom does not render the portal content, these checks pass whether or not the code shows the toast, and no test proves that a toast appears. The "Example files" dialog, also in a portal, has no test that clicks through it.

## Evidence

- Checks for the absence of the toast "An earlier upload cannot be removed" [packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L154](../../packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L154), [#L212](../../packages/app-ui/src/analysis/__tests__/test_input_upload.test.tsx#L212)
- Toasts render in `ToastPortal` [packages/app-ui/src/ui/toast.tsx#L23-L25](../../packages/app-ui/src/ui/toast.tsx#L23-L25); no test in `packages/app-ui` asserts that a toast is present

> [!IMPORTANT]
> **Investigation required.** Write a test that triggers a known toast under the same `Toaster` setup and asserts that its text is found. If the text is not found, the absence checks are vacuous.

## Fix direction

Make portal content reachable in the tests (render the portal container into the test document, or query the portal root), then pair every absence check with a test where the toast does appear, and add a click test of the "Example files" dialog.

## Validation

- Each toast absence check fails when the toast is forced to appear
