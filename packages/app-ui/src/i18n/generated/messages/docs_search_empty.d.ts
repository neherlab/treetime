export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Search_EmptyInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Search_EmptyInputs */
/**
* | output |
* | --- |
* | "No pages match the search." |
*
* @param {Docs_Search_EmptyInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_search_empty: ((inputs?: Docs_Search_EmptyInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_EmptyInputs, {
    locale?: "en" | "de";
}, {}>;
