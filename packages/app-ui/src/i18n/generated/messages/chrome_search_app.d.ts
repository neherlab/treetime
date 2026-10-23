export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_Search_AppInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Chrome_Search_AppInputs */
/**
* | output |
* | --- |
* | "Search runs, settings, examples" |
*
* @param {Chrome_Search_AppInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_search_app: ((inputs?: Chrome_Search_AppInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Search_AppInputs, {
    locale?: "en" | "de";
}, {}>;
