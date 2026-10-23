export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Site_NameInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Site_NameInputs */
/**
* | output |
* | --- |
* | "TreeTime docs" |
*
* @param {Docs_Site_NameInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_site_name: ((inputs?: Docs_Site_NameInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Site_NameInputs, {
    locale?: "en" | "de";
}, {}>;
