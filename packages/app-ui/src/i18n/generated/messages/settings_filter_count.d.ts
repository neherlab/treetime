export type LocalizedString = import('../runtime.js').LocalizedString;
export type Settings_Filter_CountInputs = {
    count: NonNullable<unknown>;
};
/**
* | countPlural | output |
* | --- | --- |
* | "one" | "{count} setting matches" |
* | "other" | "{count} settings match" |
*
* @param {Settings_Filter_CountInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export declare const settings_filter_count: ((inputs: Settings_Filter_CountInputs, options?: {
    locale?: "en" | "de";
}) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Filter_CountInputs, {
    locale?: "en" | "de";
}, {}>;
