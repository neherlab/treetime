/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Fallback_TitleInputs */



/**
* | output |
* | --- |
* | "This page is not translated yet" |
*
* @param {Docs_Fallback_TitleInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_fallback_title = /** @type {((inputs?: Docs_Fallback_TitleInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Fallback_TitleInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Diese Seite ist noch nicht übersetzt")
	return /** @type {LocalizedString} */ ("This page is not translated yet")
});