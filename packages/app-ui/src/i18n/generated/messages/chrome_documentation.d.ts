export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_DocumentationInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Chrome_DocumentationInputs */
/**
* | output |
* | --- |
* | "Documentation" |
*
* @param {Chrome_DocumentationInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_documentation: ((inputs?: Chrome_DocumentationInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_DocumentationInputs, {
    locale?: "en" | "de";
}, {}>;
