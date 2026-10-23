/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Search_UnavailableInputs */



/**
* | output |
* | --- |
* | "Search works on the built site. Run `just build-docs` and `just serve-docs` to use it." |
*
* @param {Docs_Search_UnavailableInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_search_unavailable = /** @type {((inputs?: Docs_Search_UnavailableInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Search_UnavailableInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Die Suche funktioniert auf der gebauten Seite. Führen Sie `just build-docs` und `just serve-docs` aus, um sie zu nutzen.")
	return /** @type {LocalizedString} */ ("Search works on the built site. Run `just build-docs` and `just serve-docs` to use it.")
});