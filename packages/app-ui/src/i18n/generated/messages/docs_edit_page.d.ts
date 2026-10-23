export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Edit_PageInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Edit_PageInputs */
/**
* | output |
* | --- |
* | "Edit this page on GitHub" |
*
* @param {Docs_Edit_PageInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_edit_page: ((inputs?: Docs_Edit_PageInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Edit_PageInputs, {
    locale?: "en" | "de";
}, {}>;
