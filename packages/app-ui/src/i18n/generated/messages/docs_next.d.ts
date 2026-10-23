export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_NextInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_NextInputs */
/**
* | output |
* | --- |
* | "Next" |
*
* @param {Docs_NextInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_next: ((inputs?: Docs_NextInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_NextInputs, {
    locale?: "en" | "de";
}, {}>;
