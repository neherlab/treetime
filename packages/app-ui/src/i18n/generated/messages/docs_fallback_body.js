/* eslint-disable */
import { getLocale, experimentalStaticLocale } from '../runtime.js';

/** @typedef {import('../runtime.js').LocalizedString} LocalizedString */

/** @typedef {{}} Docs_Fallback_BodyInputs */



/**
* | output |
* | --- |
* | "It is shown in English. Setting texts, help texts, and the command-line reference stay in English on every page." |
*
* @param {Docs_Fallback_BodyInputs} inputs
* @param {{ locale?: "en" | "de" }} options
* @returns {LocalizedString}
*/
export const docs_fallback_body = /** @type {((inputs?: Docs_Fallback_BodyInputs, options?: { locale?: "en" | "de" }) => LocalizedString) & import('../runtime.js').MessageMetadata<Docs_Fallback_BodyInputs, { locale?: "en" | "de" }, {}>} */ ((inputs = {}, options = {}) => {
	void inputs;
	const locale = experimentalStaticLocale ?? options.locale ?? getLocale()
	if (locale === "de") return /** @type {LocalizedString} */ ("Sie wird auf Englisch angezeigt. Einstellungstexte, Hilfetexte und die Befehlszeilenreferenz bleiben auf allen Seiten englisch.")
	return /** @type {LocalizedString} */ ("It is shown in English. Setting texts, help texts, and the command-line reference stay in English on every page.")
});