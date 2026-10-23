export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Open_NavigationInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Docs_Open_NavigationInputs */
/**
* | output |
* | --- |
* | "Open the navigation" |
*
* @param {Docs_Open_NavigationInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_open_navigation: ((inputs?: Docs_Open_NavigationInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Open_NavigationInputs, {
    locale?: "en" | "de";
}, {}>;
