export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Not_Found_BodyInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Not_Found_BodyInputs */
/**
* | output |
* | --- |
* | "No page of the docs has this address." |
*
* @param {Docs_Not_Found_BodyInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_not_found_body: ((inputs?: Docs_Not_Found_BodyInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Not_Found_BodyInputs, {
    locale?: "en" | "de";
}, {}>;
