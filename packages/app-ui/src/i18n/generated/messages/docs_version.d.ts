export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_VersionInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_VersionInputs */
/**
* | output |
* | --- |
* | "Docs version" |
*
* @param {Docs_VersionInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_version: ((inputs?: Docs_VersionInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_VersionInputs, {
    locale?: "en" | "de";
}, {}>;
