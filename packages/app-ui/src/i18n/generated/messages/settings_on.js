/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Settings_OnInputs */



/**
* | output |
* | --- |
* | "on" |
*
* @param {Settings_OnInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const settings_on = /** @type {((inputs?: Settings_OnInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Settings_OnInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("an")
	return /** @type {LocalizedString} */ ("on")
});