/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Chrome_Open_AppInputs */



/**
* | output |
* | --- |
* | "Open TreeTime" |
*
* @param {Chrome_Open_AppInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const chrome_open_app = /** @type {((inputs?: Chrome_Open_AppInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Open_AppInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("TreeTime öffnen")
	return /** @type {LocalizedString} */ ("Open TreeTime")
});