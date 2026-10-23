export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Locale_EnInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Locale_EnInputs */
/**
* | output |
* | --- |
* | "English" |
*
* @param {Docs_Locale_EnInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_locale_en: ((inputs?: Docs_Locale_EnInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Locale_EnInputs, {
    locale?: "en" | "de";
}, {}>;
