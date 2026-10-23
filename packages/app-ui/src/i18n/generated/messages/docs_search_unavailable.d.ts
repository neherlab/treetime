export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Search_UnavailableInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Search_UnavailableInputs */
/**
* | output |
* | --- |
* | "Search works on the built site. Run `just build-docs` and `just serve-docs` to use it." |
*
* @param {Docs_Search_UnavailableInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_search_unavailable: ((inputs?: Docs_Search_UnavailableInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_UnavailableInputs, {
    locale?: "en" | "de";
}, {}>;
