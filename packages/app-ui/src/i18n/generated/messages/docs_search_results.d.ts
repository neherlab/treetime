export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Search_ResultsInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Search_ResultsInputs */
/**
* | output |
* | --- |
* | "Search results" |
*
* @param {Docs_Search_ResultsInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_search_results: ((inputs?: Docs_Search_ResultsInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_ResultsInputs, {
    locale?: "en" | "de";
}, {}>;
