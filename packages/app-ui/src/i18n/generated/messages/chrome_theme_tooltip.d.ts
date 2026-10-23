export type LocalizedString = import('../runtime.js').LocalizedString;
export type Chrome_Theme_TooltipInputs = {
    label: NonNullable<unknown>;
};
/**
* | output |
* | --- |
* | "{label}. Click to change." |
*
* @param {Chrome_Theme_TooltipInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const chrome_theme_tooltip: ((inputs: Chrome_Theme_TooltipInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Theme_TooltipInputs, {
    locale?: "en" | "de";
}, {}>;
