export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Fallback_TitleInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Fallback_TitleInputs */
/**
* | output |
* | --- |
* | "This page is not translated yet" |
*
* @param {Docs_Fallback_TitleInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_fallback_title: ((inputs?: Docs_Fallback_TitleInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Fallback_TitleInputs, {
    locale?: "en" | "de";
}, {}>;
