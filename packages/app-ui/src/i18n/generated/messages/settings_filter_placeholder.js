/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Settings_Filter_PlaceholderInputs */



/**
* | output |
* | --- |
* | "Label, flag, config key, or help text" |
*
* @param {Settings_Filter_PlaceholderInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const settings_filter_placeholder = /** @type {((inputs?: Settings_Filter_PlaceholderInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_Filter_PlaceholderInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Bezeichnung, Option, Konfigurationsschlüssel oder Hilfetext")
	return /** @type {LocalizedString} */ ("Label, flag, config key, or help text")
});