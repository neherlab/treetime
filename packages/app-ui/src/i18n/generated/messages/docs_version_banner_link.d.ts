export type LocalizedString = import('../runtime.js').LocalizedString;
export type Docs_Version_Banner_LinkInputs = {
    preferred: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "Open this page in {preferred}" |
*
* @param {Docs_Version_Banner_LinkInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const docs_version_banner_link: ((inputs: Docs_Version_Banner_LinkInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Version_Banner_LinkInputs, {
    locale?: "en" | "de";
}, {}>;
