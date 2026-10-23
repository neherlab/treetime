export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Locale_DeInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Locale_DeInputs */
/**
* | output |
* | --- |
* | "Deutsch" |
*
* @param {Docs_Locale_DeInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_locale_de: ((inputs?: Docs_Locale_DeInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Locale_DeInputs, {
    locale?: "en" | "de";
}, {}>;
