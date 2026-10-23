/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Chrome_Search_AppInputs */



/**
* | output |
* | --- |
* | "Search runs, settings, examples" |
*
* @param {Chrome_Search_AppInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const chrome_search_app = /** @type {((inputs?: Chrome_Search_AppInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Chrome_Search_AppInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Läufe, Einstellungen, Beispiele suchen")
	return /** @type {LocalizedString} */ ("Search runs, settings, examples")
});