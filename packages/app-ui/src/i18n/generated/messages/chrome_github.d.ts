export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_GithubInputs = {};
/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */
/** @typedef {{}} Chrome_GithubInputs */
/**
* | output |
* | --- |
* | "TreeTime on GitHub" |
*
* @param {Chrome_GithubInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_github: ((inputs?: Chrome_GithubInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_GithubInputs, {
    locale?: "en" | "de";
}, {}>;
