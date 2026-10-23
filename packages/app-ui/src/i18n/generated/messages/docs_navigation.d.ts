export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_NavigationInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_NavigationInputs */
/**
* | output |
* | --- |
* | "Docs navigation" |
*
* @param {Docs_NavigationInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_navigation: ((inputs?: Docs_NavigationInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_NavigationInputs, {
    locale?: "en" | "de";
}, {}>;
