export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_LanguageInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_LanguageInputs */
/**
* | output |
* | --- |
* | "Language" |
*
* @param {Docs_LanguageInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_language: ((inputs?: Docs_LanguageInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_LanguageInputs, {
    locale?: "en" | "de";
}, {}>;
