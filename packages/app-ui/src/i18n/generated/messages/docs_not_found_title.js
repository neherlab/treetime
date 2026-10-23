/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Not_Found_TitleInputs */



/**
* | output |
* | --- |
* | "Page not found" |
*
* @param {Docs_Not_Found_TitleInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_not_found_title = /** @type {((inputs?: Docs_Not_Found_TitleInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Not_Found_TitleInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Seite nicht gefunden")
	return /** @type {LocalizedString} */ ("Page not found")
});