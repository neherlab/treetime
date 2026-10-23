/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Settings_Filter_LabelInputs */



/**
* | output |
* | --- |
* | "Filter settings" |
*
* @param {Settings_Filter_LabelInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const settings_filter_label = /** @type {((inputs?: Settings_Filter_LabelInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Filter_LabelInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Einstellungen filtern")
	return /** @type {LocalizedString} */ ("Filter settings")
});