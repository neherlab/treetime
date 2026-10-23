export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_PreviousInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_PreviousInputs */
/**
* | output |
* | --- |
* | "Previous" |
*
* @param {Docs_PreviousInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_previous: ((inputs?: Docs_PreviousInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_PreviousInputs, {
    locale?: "en" | "de";
}, {}>;
