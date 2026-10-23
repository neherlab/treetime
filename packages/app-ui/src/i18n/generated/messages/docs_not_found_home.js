/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Not_Found_HomeInputs */



/**
* | output |
* | --- |
* | "Go to the docs home page" |
*
* @param {Docs_Not_Found_HomeInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_not_found_home = /** @type {((inputs?: Docs_Not_Found_HomeInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Not_Found_HomeInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Zur Startseite der Doku")
	return /** @type {LocalizedString} */ ("Go to the docs home page")
});