export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Not_Found_HomeInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Not_Found_HomeInputs */
/**
* | output |
* | --- |
* | "Go to the docs home page" |
*
* @param {Docs_Not_Found_HomeInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_not_found_home: ((inputs?: Docs_Not_Found_HomeInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Not_Found_HomeInputs, {
    locale?: "en" | "de";
}, {}>;
