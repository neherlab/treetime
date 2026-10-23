/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Chrome_Theme_DarkInputs */



/**
* | output |
* | --- |
* | "Dark theme" |
*
* @param {Chrome_Theme_DarkInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const chrome_theme_dark = /** @type {((inputs?: Chrome_Theme_DarkInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Theme_DarkInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Dunkles Design")
	return /** @type {LocalizedString} */ ("Dark theme")
});