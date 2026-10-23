export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Search_LoadingInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Search_LoadingInputs */
/**
* | output |
* | --- |
* | "Loading the search index" |
*
* @param {Docs_Search_LoadingInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_search_loading: ((inputs?: Docs_Search_LoadingInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_LoadingInputs, {
    locale?: "en" | "de";
}, {}>;
