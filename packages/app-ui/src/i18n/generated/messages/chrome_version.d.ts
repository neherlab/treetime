export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_VersionInputs = {
    version: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "v{version}" |
*
* @param {Chrome_VersionInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_version: ((inputs: Chrome_VersionInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_VersionInputs, {
    locale?: "en" | "de";
}, {}>;
