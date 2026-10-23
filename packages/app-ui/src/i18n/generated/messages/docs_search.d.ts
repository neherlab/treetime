export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_SearchInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_SearchInputs */
/**
* | output |
* | --- |
* | "Search the docs" |
*
* @param {Docs_SearchInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_search: ((inputs?: Docs_SearchInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_SearchInputs, {
    locale?: "en" | "de";
}, {}>;
