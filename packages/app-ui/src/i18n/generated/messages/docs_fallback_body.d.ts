export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Fallback_BodyInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Fallback_BodyInputs */
/**
* | output |
* | --- |
* | "It is shown in English. Setting texts, help texts, and the command-line reference stay in English on every page." |
*
* @param {Docs_Fallback_BodyInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_fallback_body: ((inputs?: Docs_Fallback_BodyInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Fallback_BodyInputs, {
    locale?: "en" | "de";
}, {}>;
