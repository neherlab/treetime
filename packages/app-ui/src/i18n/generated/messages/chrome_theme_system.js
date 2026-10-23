/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Chrome_Theme_SystemInputs */



/**
* | output |
* | --- |
* | "Theme follows the system" |
*
* @param {Chrome_Theme_SystemInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const chrome_theme_system = /** @type {((inputs?: Chrome_Theme_SystemInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Theme_SystemInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Design folgt dem System")
	return /** @type {LocalizedString} */ ("Theme follows the system")
});