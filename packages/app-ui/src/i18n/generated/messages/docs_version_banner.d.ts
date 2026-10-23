export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Version_BannerInputs = {
    current: NonNullable<unknown>;
    preferred: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "These are the docs of TreeTime {current}. The recommended version is {preferred}." |
*
* @param {Docs_Version_BannerInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_version_banner: ((inputs: Docs_Version_BannerInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Version_BannerInputs, {
    locale?: "en" | "de";
}, {}>;
