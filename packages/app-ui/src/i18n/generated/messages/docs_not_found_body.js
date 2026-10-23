/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Not_Found_BodyInputs */



/**
* | output |
* | --- |
* | "No page of the docs has this address." |
*
* @param {Docs_Not_Found_BodyInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_not_found_body = /** @type {((inputs?: Docs_Not_Found_BodyInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Not_Found_BodyInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Keine Seite der Doku hat diese Adresse.")
	return /** @type {LocalizedString} */ ("No page of the docs has this address.")
});