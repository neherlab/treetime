export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Search_PlaceholderInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Search_PlaceholderInputs */
/**
* | output |
* | --- |
* | "Search pages and headings" |
*
* @param {Docs_Search_PlaceholderInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_search_placeholder: ((inputs?: Docs_Search_PlaceholderInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_PlaceholderInputs, {
    locale?: "en" | "de";
}, {}>;
