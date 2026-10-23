export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_TocInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_TocInputs */
/**
* | output |
* | --- |
* | "On this page" |
*
* @param {Docs_TocInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_toc: ((inputs?: Docs_TocInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_TocInputs, {
    locale?: "en" | "de";
}, {}>;
