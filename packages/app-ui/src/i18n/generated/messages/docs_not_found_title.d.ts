export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Not_Found_TitleInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Not_Found_TitleInputs */
/**
* | output |
* | --- |
* | "Page not found" |
*
* @param {Docs_Not_Found_TitleInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_not_found_title: ((inputs?: Docs_Not_Found_TitleInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Not_Found_TitleInputs, {
    locale?: "en" | "de";
}, {}>;
