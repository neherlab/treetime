/**
 * The project's base locale.
 *
 * @example
 *   if (locale === baseLocale) {
 *     // do something
 *   }
 */
export declare const baseLocale = "en";
/**
 * The project's locales that have been specified in the settings.
 *
 * @example
 *   if (locales.includes(userSelectedLocale) === false) {
 *     throw new Error('Locale is not available');
 *   }
 */
export declare const locales: readonly ["en", "de"];
/** @type {string} */
export declare const cookieName: string;
/** @type {number} */
export declare const cookieMaxAge: number;
/** @type {string} */
export declare const cookieDomain: string;
export declare const cookieSameSite: "lax" | "strict" | "none" | undefined;
export declare const cookieSecure: boolean | "auto";
/** @type {string} */
export declare const localStorageKey: string;
/**
 * @type {Array<"cookie" | "baseLocale" | "globalVariable" | "url" | "preferredLanguage" | "localStorage" | `custom-${string}`>}
 */
export declare const strategy: Array<"cookie" | "baseLocale" | "globalVariable" | "url" | "preferredLanguage" | "localStorage" | `custom-${string}`>;
/**
 * Route-level strategy overrides.
 *
 * `match` uses URLPattern syntax.
 *
 * @type {Array<{
 *   match: string;
 *   strategy?: Array<"cookie" | "baseLocale" | "globalVariable" | "url" | "preferredLanguage" | "localStorage" | `custom-${string}`>;
 *   exclude?: boolean;
 * }>}
 */
export declare const routeStrategies: Array<{
    match: string;
    strategy?: Array<"cookie" | "baseLocale" | "globalVariable" | "url" | "preferredLanguage" | "localStorage" | `custom-${string}`>;
    exclude?: boolean;
}>;
/**
 * The used URL patterns.
 *
 * @type {Array<{ pattern: string, localized: Array<[Locale, string]> }>}
 */
export declare const urlPatterns: Array<{
    pattern: string;
    localized: Array<[Locale, string]>;
}>;
/**
 * Controls trailing slash canonicalization for localized URLs.
 *
 * @type {"always" | "never" | undefined}
 */
export declare const trailingSlash: "always" | "never" | undefined;
export type ParaglideAsyncLocalStorage = {
    getStore(): {
        locale?: Locale;
        origin?: string;
        messageCalls?: Set<string>;
    } | undefined;
    run: (store: {
        locale?: Locale;
        origin?: string;
        messageCalls?: Set<string>;
    }, cb: any) => any;
};
/**
 * @typedef {{
 * 		getStore(): {
 *   		locale?: Locale,
 * 			origin?: string,
 * 			messageCalls?: Set<string>
 *   	} | undefined,
 * 		run: (store: { locale?: Locale, origin?: string, messageCalls?: Set<string>},
 *    cb: any) => any
 * }} ParaglideAsyncLocalStorage
 */
/**
 * Server side async local storage that is set by `serverMiddleware()`.
 *
 * The variable is used to retrieve the locale and origin in a server-side
 * rendering context without effecting other requests.
 *
 * @type {ParaglideAsyncLocalStorage | undefined}
 */
export declare let serverAsyncLocalStorage: ParaglideAsyncLocalStorage | undefined;
/**
 * Returns the current server-side async local storage instance.
 *
 * Accessing the mutable value through a function keeps it observable when
 * module interceptors wrap exported bindings and snapshot their initial value.
 *
 * @returns {ParaglideAsyncLocalStorage | undefined}
 */
export declare function getServerAsyncLocalStorage(): ParaglideAsyncLocalStorage | undefined;
export declare const disableAsyncLocalStorage = false;
export declare const experimentalMiddlewareLocaleSplitting = false;
export declare const isServer: boolean;
/** @type {Locale | undefined} */
export declare const experimentalStaticLocale: Locale | undefined;
/**
 * Sets the server side async local storage.
 *
 * The function is needed because the `runtime.js` file
 * must define the `serverAsyncLocalStorage` variable to
 * avoid a circular import between `runtime.js` and
 * `server.js` files.
 *
 * @param {ParaglideAsyncLocalStorage | undefined} value
 */
export declare function overwriteServerAsyncLocalStorage(value: ParaglideAsyncLocalStorage | undefined): void;
/**
 * Get the current locale.
 *
 * The locale is resolved using your configured strategies (URL, cookie, localStorage, etc.)
 * in the order they are defined. In SSR contexts, the locale is retrieved from AsyncLocalStorage
 * which is set by the `paraglideMiddleware()`.
 *
 * @see https://paraglidejs.com/strategy - Configure locale detection strategies
 *
 * @example
 *   if (getLocale() === 'de') {
 *     console.log('Germany 🇩🇪');
 *   } else if (getLocale() === 'nl') {
 *     console.log('Netherlands 🇳🇱');
 *   }
 *
 * @returns {Locale} The current locale.
 */
export declare let getLocale: () => Locale;
/**
 * Resolve locale for a given URL using route-aware strategies.
 *
 * @param {string | URL} url
 * @returns {Locale}
 */
export declare function getLocaleForUrl(url: string | URL): Locale;
/**
 * Overwrite the `getLocale()` function.
 *
 * Use this function to overwrite how the locale is resolved. This is useful
 * for custom locale resolution or advanced use cases like SSG with concurrent rendering.
 *
 * @see https://paraglidejs.com/strategy
 *
 * @example
 *   overwriteGetLocale(() => {
 *     return Cookies.get('locale') ?? baseLocale
 *   });
 *
 * @param {() => Locale} fn - The new implementation for `getLocale()`.
 */
export declare const overwriteGetLocale: (fn: () => Locale) => void;
/**
 * Get writing direction for a locale.
 *
 * Uses `Intl.Locale` text info when available and falls back to a
 * language-based RTL check for runtimes without `getTextInfo()`.
 *
 * @example
 *   getTextDirection(); // "ltr" or "rtl" for current locale
 *   getTextDirection("ar"); // "rtl"
 *   getTextDirection("en"); // "ltr"
 *
 * @param {string} [locale] - Target locale. If not provided, uses `getLocale()`
 * @returns {"ltr" | "rtl"}
 */
export declare function getTextDirection(locale?: string): "ltr" | "rtl";
export type SetLocaleFn = (newLocale: Locale, options?: {
    reload?: boolean;
}) => void | Promise<void>;
/**
 * @typedef {(newLocale: Locale, options?: { reload?: boolean }) => void | Promise<void>} SetLocaleFn
 */
/**
 * Set the locale.
 *
 * Updates the locale using your configured strategies (cookie, localStorage, URL, etc.).
 * By default, this navigates the client to the localized URL or reloads the current
 * document to reflect the new locale. `reload: false` is a narrow browser-only escape
 * hatch for a fully client-rendered, non-URL-routed surface that owns its reactive
 * updates and document state. It does not re-render the UI or update the document.
 * Do not use it for normal locale pickers, URL-routed pages, or switching an SSR,
 * SSG, or hydrated document. It is incompatible with per-locale builds.
 *
 * If any custom strategy's `setLocale` function is async, then this function
 * will become async as well.
 *
 * @see https://paraglidejs.com/strategy
 *
 * @example
 *   setLocale('en');
 *
 * @example
 *   setLocale('en', { reload: false });
 *
 * @type {SetLocaleFn}
 */
export declare let setLocale: SetLocaleFn;
/**
 * Overwrite the `setLocale()` function.
 *
 * Use this function to overwrite how the locale is set. For example,
 * modify a cookie, env variable, or a user's preference.
 *
 * @example
 *   overwriteSetLocale((newLocale) => {
 *     // set the locale in a cookie
 *     return Cookies.set('locale', newLocale)
 *   });
 *
 * @param {SetLocaleFn} fn
 */
export declare const overwriteSetLocale: (fn: SetLocaleFn) => void;
/**
 * The origin of the current URL.
 *
 * Defaults to "http://example.com" in non-browser environments. If this
 * behavior is not desired, the implementation can be overwritten
 * by `overwriteGetUrlOrigin()`.
 *
 * @type {() => string}
 */
export declare let getUrlOrigin: () => string;
/**
 * Overwrite the getUrlOrigin function.
 *
 * Use this function in server environments to
 * define how the URL origin is resolved.
 *
 * @param {() => string} fn - The new implementation for `getUrlOrigin()`.
 */
export declare let overwriteGetUrlOrigin: (fn: () => string) => void;
/**
 * Coerces a locale-like string to the canonical locale value used by the runtime.
 *
 * @param {unknown} value
 * @returns {Locale | undefined}
 */
export declare function toLocale(value: unknown): Locale | undefined;
/**
 * Check if something is an available locale with the canonical project casing.
 *
 * @example
 *   if (isLocale(params.locale)) {
 *     setLocale(params.locale);
 *   } else {
 *     setLocale('en');
 *   }
 *
 * Use `toLocale()` when you want case-insensitive matching and canonicalization.
 *
 * @param {unknown} locale
 * @returns {locale is Locale}
 */
export declare function isLocale(locale: unknown): locale is Locale;
/**
 * Asserts that the input can be normalized to a locale.
 *
 * @param {unknown} input - The input to check.
 * @returns {Locale} The input normalized to a Locale.
 * @throws {Error} If the input is not a locale.
 */
export declare function assertIsLocale(input: unknown): Locale;
export type ExtractLocaleFromRequestOptions = {
    /**
     * - Effective request URL to use for route matching and locale detection with the URL strategy.
     */
    effectiveRequestUrl?: string | URL;
};
/**
 * @typedef {object} ExtractLocaleFromRequestOptions
 * @property {string | URL} [effectiveRequestUrl] - Effective request URL to use for route matching and locale detection with the URL strategy.
 */
/**
 * Extracts a locale from a request.
 *
 * Use the function on the server to extract the locale
 * from a request.
 *
 * The function goes through the strategies in the order
 * they are defined. If a strategy returns an invalid locale,
 * it will fall back to the next strategy.
 *
 * Note: Custom server strategies are not supported in this synchronous version.
 * Use `extractLocaleFromRequestAsync` if you need custom server strategies with async getLocale methods.
 *
 * @example
 *   const locale = extractLocaleFromRequest(request);
 *
 * @param {Request} request
 * @param {ExtractLocaleFromRequestOptions} [options]
 * @returns {Locale}
 */
export declare const extractLocaleFromRequest: (request: Request, options?: ExtractLocaleFromRequestOptions) => Locale;
/**
 * Extracts a locale from a request using the provided strategy order.
 *
 * @param {Request} request
 * @param {typeof strategy} strategies
 * @param {string | URL} [url]
 * @returns {Locale}
 */
export declare const extractLocaleFromRequestWithStrategies: (request: Request, strategies: typeof strategy, url?: string | URL) => Locale;
/**
 * Asynchronously extracts a locale from a request.
 *
 * This function supports async custom server strategies, unlike the synchronous
 * `extractLocaleFromRequest`. Use this function when you have custom server strategies
 * that need to perform asynchronous operations (like database calls) in their getLocale method.
 *
 * The function first processes any custom server strategies asynchronously, then falls back
 * to the synchronous `extractLocaleFromRequest` for all other strategies.
 *
 * @see {@link https://github.com/opral/inlang-paraglide-js/issues/527#issuecomment-2978151022}
 *
 * @example
 *   // Basic usage
 *   const locale = await extractLocaleFromRequestAsync(request);
 *
 * @example
 *   // With custom async server strategy
 *   defineCustomServerStrategy("custom-database", {
 *     getLocale: async (request) => {
 *       const userId = extractUserIdFromRequest(request);
 *       return await getUserLocaleFromDatabase(userId);
 *     }
 *   });
 *
 *   const locale = await extractLocaleFromRequestAsync(request);
 *
 * @param {Request} request - The request object to extract the locale from.
 * @param {{ effectiveRequestUrl?: string | URL }} [options] - Effective request URL to use for route matching and locale detection with the URL strategy.
 * @returns {Promise<Locale>} The extracted locale.
 */
export declare const extractLocaleFromRequestAsync: (request: Request, options?: {
    effectiveRequestUrl?: string | URL;
}) => Promise<Locale>;
/**
 * Extracts a cookie from the document.
 *
 * Will return undefined if the document is not available or if the cookie is not set.
 * The `document` object is not available in server-side rendering, so this function should not be called in that context.
 *
 * @returns {Locale | undefined}
 */
export declare function extractLocaleFromCookie(): Locale | undefined;
/**
 * Extracts a locale from the accept-language header.
 *
 * Use the function on the server to extract the locale
 * from the accept-language header that is sent by the client.
 *
 * @example
 *   const locale = extractLocaleFromHeader(request);
 *
 * @param {Request} request - The request object to extract the locale from.
 * @returns {Locale | undefined} The negotiated preferred language.
 */
export declare function extractLocaleFromHeader(request: Request): Locale | undefined;
/**
 * Negotiates a preferred language from navigator.languages.
 *
 * Use the function on the client to extract the locale
 * from the navigator.languages array.
 *
 * @example
 *   const locale = extractLocaleFromNavigator();
 *
 * @returns {Locale | undefined}
 */
export declare function extractLocaleFromNavigator(): Locale | undefined;
/**
 * Extracts the locale from a given URL using native URLPattern.
 *
 * The built-in default `/:locale/...` routing is case-insensitive because it
 * canonicalizes the first path segment with `toLocale()`. Custom `urlPatterns`
 * keep URLPattern's normal exact matching semantics for path segments.
 *
 * @param {URL|string} url - The full URL from which to extract the locale.
 * @returns {Locale|undefined} The extracted locale, or undefined if no locale is found.
 */
export declare function extractLocaleFromUrl(url: URL | string): Locale | undefined;
/**
 * Lower-level URL localization function, primarily used in server contexts.
 *
 * This function is designed for server-side usage where you need precise control
 * over URL localization, such as in middleware or request handlers. It works with
 * URL objects and always returns absolute URLs.
 *
 * For client-side UI components, use `localizeHref()` instead, which provides
 * a more convenient API with relative paths and automatic locale detection.
 *
 * @see https://paraglidejs.com/i18n-routing
 *
 * @example
 * ```typescript
 * // Server middleware example
 * app.use((req, res, next) => {
 *   const url = new URL(req.url, `${req.protocol}://${req.headers.host}`);
 *   const localized = localizeUrl(url, { locale: "de" });
 *
 *   if (localized.href !== url.href) {
 *     return res.redirect(localized.href);
 *   }
 *   next();
 * });
 * ```
 *
 * @example
 * ```typescript
 * // Using with URL patterns
 * const url = new URL("https://example.com/about");
 * localizeUrl(url, { locale: "de" });
 * // => URL("https://example.com/de/about")
 *
 * // Using with domain-based localization
 * const url = new URL("https://example.com/store");
 * localizeUrl(url, { locale: "de" });
 * // => URL("https://de.example.com/store")
 * ```
 *
 * @param {string | { readonly href: string }} url - The URL to localize. Accepts readonly URL objects. If string, must be absolute.
 * @param {object} [options] - Options for localization
 * @param {Locale} [options.locale] - Target locale. If not provided, uses getLocale()
 * @returns {URL} The localized URL, always absolute
 */
export declare function localizeUrl(url: string | {
    readonly href: string;
}, options?: {
    locale?: Locale;
}): URL;
/**
 * Low-level URL de-localization function, primarily used in server contexts.
 *
 * This function is designed for server-side usage where you need precise control
 * over URL de-localization, such as in middleware or request handlers. It works with
 * URL objects and always returns absolute URLs.
 *
 * For client-side UI components, use `deLocalizeHref()` instead, which provides
 * a more convenient API with relative paths.
 *
 * @see https://paraglidejs.com/i18n-routing
 *
 * @example
 * ```typescript
 * // Server middleware example
 * app.use((req, res, next) => {
 *   const url = new URL(req.url, `${req.protocol}://${req.headers.host}`);
 *   const baseUrl = deLocalizeUrl(url);
 *
 *   // Store the base URL for later use
 *   req.baseUrl = baseUrl;
 *   next();
 * });
 * ```
 *
 * @example
 * ```typescript
 * // Using with URL patterns
 * const url = new URL("https://example.com/de/about");
 * deLocalizeUrl(url); // => URL("https://example.com/about")
 *
 * // Using with domain-based localization
 * const url = new URL("https://de.example.com/store");
 * deLocalizeUrl(url); // => URL("https://example.com/store")
 * ```
 *
 * @param {string | { readonly href: string }} url - The URL to de-localize. Accepts readonly URL objects. If string, must be absolute.
 * @returns {URL} The de-localized URL, always absolute
 */
export declare function deLocalizeUrl(url: string | {
    readonly href: string;
}): URL;
/**
 * Aggregates named groups from various parts of the URLPattern match result.
 *
 *
 * @param {any} match - The URLPattern match result object.
 * @returns {Record<string, string | null | undefined>} An object containing all named groups from the match.
 */
export declare function aggregateGroups(match: any): Record<string, string | null | undefined>;
export type FastPathPattern = {
    protocol: string | undefined;
    hostname: string | undefined;
    port: string | undefined;
    pathnamePrefix: string;
    pathMode: "segments" | "catch-all-optional" | "catch-all-required";
};
export type FastPathRoute = {
    base: FastPathPattern;
    localized: Array<{
        locale: string;
        pattern: FastPathPattern;
    }>;
};
/**
 * Match route policy against both the public URL and its canonical URL.
 *
 * The function is deliberately separate from variables.js: configuration is
 * inert data, while canonicalization and route selection form a routing layer.
 *
 * @param {string | URL} url
 * @returns {{ match: string; strategy?: typeof strategy; exclude?: boolean } | undefined}
 */
export declare function findMatchingRouteStrategy(url: string | URL): {
    match: string;
    strategy?: typeof strategy;
    exclude?: boolean;
} | undefined;
/**
 * Returns the strategy to use for a specific URL.
 *
 * If no route strategy matches (or the matching rule is `exclude: true`),
 * the global strategy is returned.
 *
 * @param {string | URL} url
 * @returns {typeof strategy}
 */
export declare function getStrategyForUrl(url: string | URL): typeof strategy;
/**
 * Returns whether the given URL is excluded from middleware i18n processing.
 *
 * @param {string | URL} url
 * @returns {boolean}
 */
export declare function isExcludedByRouteStrategy(url: string | URL): boolean;
export type ShouldRedirectServerInput = {
    request: Request;
    /**
     * - Effective request URL to use for route matching, locale detection with the URL strategy, and redirect targets.
     */
    effectiveRequestUrl?: string | URL;
    locale?: Locale;
};
export type ShouldRedirectClientInput = {
    request?: undefined;
    url?: string | URL;
    locale?: Locale;
};
export type ShouldRedirectInput = ShouldRedirectServerInput | ShouldRedirectClientInput;
export type ShouldRedirectResult = {
    /**
     * - Indicates whether the consumer should perform a redirect.
     */
    shouldRedirect: boolean;
    /**
     * - Locale resolved using the configured strategies.
     */
    locale: Locale;
    /**
     * - Destination URL when a redirect is required.
     */
    redirectUrl: URL | undefined;
};
/**
 * @typedef {object} ShouldRedirectServerInput
 * @property {Request} request
 * @property {string | URL} [effectiveRequestUrl] - Effective request URL to use for route matching, locale detection with the URL strategy, and redirect targets.
 * @property {Locale} [locale]
 *
 * @typedef {object} ShouldRedirectClientInput
 * @property {undefined} [request]
 * @property {string | URL} [url]
 * @property {Locale} [locale]
 *
 * @typedef {ShouldRedirectServerInput | ShouldRedirectClientInput} ShouldRedirectInput
 *
 * @typedef {object} ShouldRedirectResult
 * @property {boolean} shouldRedirect - Indicates whether the consumer should perform a redirect.
 * @property {Locale} locale - Locale resolved using the configured strategies.
 * @property {URL | undefined} redirectUrl - Destination URL when a redirect is required.
 */
/**
 * Determines whether a redirect is required to align the current URL with the active locale.
 *
 * This helper mirrors the logic that powers `paraglideMiddleware`, but works in both server
 * and client environments. It evaluates the configured strategies in order, computes the
 * canonical localized URL, and reports when the current URL does not match.
 *
 * When called in the browser without arguments, the current `window.location.href` is used.
 *
 * @see https://paraglidejs.com/i18n-routing#redirects
 *
 * @example
 * // Client side usage (e.g. TanStack Router beforeLoad hook)
 * async function beforeLoad({ location }) {
 *   const decision = await shouldRedirect({ url: location.href });
 *
 *   if (decision.shouldRedirect) {
 *     throw redirect({ to: decision.redirectUrl.href });
 *   }
 * }
 *
 * @example
 * // Server side usage with a Request
 * export async function handle(request) {
 *   const decision = await shouldRedirect({ request });
 *
 *   if (decision.shouldRedirect) {
 *     return Response.redirect(decision.redirectUrl, 307);
 *   }
 *
 *   return render(request, decision.locale);
 * }
 *
 * @example
 * // Server side usage behind a proxy where request.url is not public-facing
 * export async function handle(request) {
 *   const effectiveRequestUrl = new URL(request.url);
 *   effectiveRequestUrl.protocol = "https:";
 *   effectiveRequestUrl.host = "example.com";
 *
 *   const decision = await shouldRedirect({
 *     request,
 *     effectiveRequestUrl,
 *   });
 *
 *   if (decision.shouldRedirect) {
 *     return Response.redirect(decision.redirectUrl, 307);
 *   }
 * }
 *
 * @param {ShouldRedirectInput} [input]
 * @returns {Promise<ShouldRedirectResult>}
 */
export declare function shouldRedirect(input?: ShouldRedirectInput): Promise<ShouldRedirectResult>;
/**
 * High-level URL localization function optimized for client-side UI usage.
 *
 * This is a convenience wrapper around `localizeUrl()` that provides features
 * needed in UI:
 *
 * - Accepts relative paths (e.g., "/about")
 * - Returns relative paths when possible
 * - Automatically detects current locale if not specified
 * - Handles string input/output instead of URL objects
 *
 * @see https://paraglidejs.com/i18n-routing
 *
 * @example
 * ```typescript
 * // In a React/Vue/Svelte component
 * const NavLink = ({ href }) => {
 *   // Automatically uses current locale, keeps path relative
 *   return <a href={localizeHref(href)}>...</a>;
 * };
 *
 * // Examples:
 * localizeHref("/about")
 * // => "/de/about" (if current locale is "de")
 * localizeHref("/store", { locale: "fr" })
 * // => "/fr/store" (explicit locale)
 *
 * // Cross-origin links remain absolute
 * localizeHref("https://other-site.com/about")
 * // => "https://other-site.com/de/about"
 *
 * // Protocol-relative links stay protocol-relative
 * localizeHref("//other-site.com/about")
 * // => "//other-site.com/de/about"
 * ```
 *
 * For server-side URL localization (e.g., in middleware), use `localizeUrl()`
 * which provides more precise control over URL handling.
 *
 * @param {string} href - The href to localize (can be relative or absolute)
 * @param {object} [options] - Options for localization
 * @param {Locale} [options.locale] - Target locale. If not provided, uses `getLocale()`
 * @returns {string} The localized href, relative if input was relative
 */
export declare function localizeHref(href: string, options?: {
    locale?: Locale;
}): string;
/**
 * High-level URL de-localization function optimized for client-side UI usage.
 *
 * This is a convenience wrapper around `deLocalizeUrl()` that provides features
 * needed in the UI:
 *
 * - Accepts relative paths (e.g., "/de/about")
 * - Returns relative paths when possible
 * - Handles string input/output instead of URL objects
 *
 * @see https://paraglidejs.com/i18n-routing
 *
 * @example
 * ```typescript
 * // In a React/Vue/Svelte component
 * const LocaleSwitcher = ({ href }) => {
 *   // Remove locale prefix before switching
 *   const baseHref = deLocalizeHref(href);
 *   return locales.map(locale =>
 *     <a href={localizeHref(baseHref, { locale })}>
 *       Switch to {locale}
 *     </a>
 *   );
 * };
 *
 * // Examples:
 * deLocalizeHref("/de/about")  // => "/about"
 * deLocalizeHref("/fr/store")  // => "/store"
 *
 * // Cross-origin links remain absolute
 * deLocalizeHref("https://example.com/de/about")
 * // => "https://example.com/about"
 *
 * // Protocol-relative links stay protocol-relative
 * deLocalizeHref("//example.com/de/about")
 * // => "//example.com/about"
 * ```
 *
 * For server-side URL de-localization (e.g., in middleware), use `deLocalizeUrl()`
 * which provides more precise control over URL handling.
 *
 * @param {string} href - The href to de-localize (can be relative or absolute)
 * @returns {string} The de-localized href, relative if input was relative
 */
export declare function deLocalizeHref(href: string): string;
/**
 * @param {string} safeModuleId
 * @param {Locale} locale
 */
export declare function trackMessageCall(safeModuleId: string, locale: Locale): void;
/**
 * Generates localized URL variants for all provided URLs based on your configured locales and URL patterns.
 *
 * This function is essential for Static Site Generation (SSG) where you need to tell your framework
 * which pages to pre-render at build time. It's also useful for generating sitemaps and
 * `<link rel="alternate" hreflang>` tags for SEO.
 *
 * The function respects your `urlPatterns` configuration - if you have translated pathnames
 * (e.g., `/about` → `/ueber-uns` for German), it will generate the correct localized paths.
 *
 * @see https://paraglidejs.com/static-site-generation
 *
 * @example
 * // Basic usage - generate all locale variants for a list of paths
 * const localizedUrls = generateStaticLocalizedUrls([
 *   "/",
 *   "/about",
 *   "/blog/post-1",
 * ]);
 * // Returns URL objects for each locale:
 * // ["/en/", "/de/", "/en/about", "/de/about", "/en/blog/post-1", "/de/blog/post-1"]
 *
 * @example
 * // Use with framework SSG APIs
 * // SvelteKit
 * export function entries() {
 *   const paths = ["/", "/about", "/contact"];
 *   return generateStaticLocalizedUrls(paths).map(url => ({
 *     locale: extractLocaleFromUrl(url)
 *   }));
 * }
 *
 * @example
 * // Sitemap generation
 * const allPages = ["/", "/about", "/blog"];
 * const sitemapUrls = generateStaticLocalizedUrls(allPages);
 *
 * @param {(string | URL)[]} urls - List of canonical URLs or paths to generate localized versions for.
 *   Can be absolute URLs (`https://example.com/about`) or paths (`/about`).
 *   Paths are resolved against `http://localhost` internally.
 * @returns {URL[]} Array of URL objects representing all localized variants.
 *   The order follows each input URL with all its locale variants before moving to the next URL.
 */
export declare function generateStaticLocalizedUrls(urls: (string | URL)[]): URL[];
export type BuiltInStrategy = "cookie" | "baseLocale" | "globalVariable" | "url" | "preferredLanguage" | "localStorage";
export type CustomStrategy = `custom_${string}`;
export type Strategy = BuiltInStrategy | CustomStrategy;
export type Strategies = Array<Strategy>;
export type CustomServerStrategyHandler = {
    getLocale: (request?: Request) => Promise<string | undefined> | (string | undefined);
};
export type CustomClientStrategyHandler = {
    getLocale: () => Promise<string | undefined> | (string | undefined);
    setLocale: (locale: string) => Promise<void> | void;
};
/**
 * @typedef {"cookie" | "baseLocale" | "globalVariable" | "url" | "preferredLanguage" | "localStorage"} BuiltInStrategy
 */
/**
 * @typedef {`custom_${string}`} CustomStrategy
 */
/**
 * @typedef {BuiltInStrategy | CustomStrategy} Strategy
 */
/**
 * @typedef {Array<Strategy>} Strategies
 */
/**
 * @typedef {{ getLocale: (request?: Request) => Promise<string | undefined> | (string | undefined) }} CustomServerStrategyHandler
 */
/**
 * @typedef {{ getLocale: () => Promise<string|undefined> | (string | undefined), setLocale: (locale: string) => Promise<void> | void }} CustomClientStrategyHandler
 */
/** @type {Map<string, CustomServerStrategyHandler>} */
export declare const customServerStrategies: Map<string, CustomServerStrategyHandler>;
/** @type {Map<string, CustomClientStrategyHandler>} */
export declare const customClientStrategies: Map<string, CustomClientStrategyHandler>;
/**
 * Checks if the given strategy is a custom strategy.
 *
 * @param {unknown} strategy The name of the custom strategy to validate.
 * Must be a string that starts with "custom-" followed by alphanumeric characters, hyphens, or underscores.
 * @returns {boolean} Returns true if it is a custom strategy, false otherwise.
 */
export declare function isCustomStrategy(strategy: unknown): boolean;
/**
 * Defines a custom strategy that is executed on the server.
 *
 * @see https://paraglidejs.com/strategy#write-your-own-strategy
 *
 * @param {string} strategy The name of the custom strategy to define. Must follow the pattern custom-name with alphanumeric characters, hyphens, or underscores.
 * @param {CustomServerStrategyHandler} handler The handler for the custom strategy, which should implement
 * the method getLocale.
 * @returns {void}
 */
export declare function defineCustomServerStrategy(strategy: string, handler: CustomServerStrategyHandler): void;
/**
 * Defines a custom strategy that is executed on the client.
 *
 * @see https://paraglidejs.com/strategy#write-your-own-strategy
 *
 * @param {string} strategy The name of the custom strategy to define. Must follow the pattern custom-name with alphanumeric characters, hyphens, or underscores.
 * @param {CustomClientStrategyHandler} handler The handler for the custom strategy, which should implement the
 * methods getLocale and setLocale.
 * @returns {void}
 */
export declare function defineCustomClientStrategy(strategy: string, handler: CustomClientStrategyHandler): void;
export {};
export type Locale = typeof locales[number];
export type LocalizedString = string & {
    readonly __brand: 'LocalizedString';
};
export type MessageMarkupOption = {
    name: string;
    value: unknown;
};
export type MessageMarkupAttribute = {
    name: string;
    value: string | true;
};
export type MessageMarkupOptions = Record<string, unknown>;
export type MessageMarkupAttributes = Record<string, string | true>;
export type MessageMarkupTag = {
    options: MessageMarkupOptions;
    attributes: MessageMarkupAttributes;
    children: boolean;
};
export type MessageMarkupSchema = Record<string, MessageMarkupTag>;
export type MessageMetadata<Inputs, Options, Markup extends MessageMarkupSchema = MessageMarkupSchema> = {
    readonly __paraglide?: {
        inputs: Inputs;
        options: Options;
        markup: Markup;
    };
};
export type MessagePart = {
    type: "text";
    value: string;
} | {
    type: "markup-start";
    name: string;
    options: MessageMarkupOptions;
    attributes: MessageMarkupAttributes;
} | {
    type: "markup-end";
    name: string;
    options: MessageMarkupOptions;
    attributes: MessageMarkupAttributes;
} | {
    type: "markup-standalone";
    name: string;
    options: MessageMarkupOptions;
    attributes: MessageMarkupAttributes;
};
export type MessageFunction = (inputs?: Record<string, never>) => LocalizedString;
export type MessageBundleFunction<T extends string> = (params: Record<string, never>, options: {
    locale: T;
}) => LocalizedString;
/**
 * A locale that is available in the project.
 *
 * @example
 *   setLocale(request.locale as Locale)
 *
 * @typedef {typeof locales[number]} Locale
 */
/**
 * A branded type representing a localized string.
 *
 * Message functions return this type instead of \`string\`, enabling TypeScript
 * to distinguish translated strings from regular strings at compile time.
 * This allows you to enforce that only properly localized content is used
 * in your UI components.
 *
 * Since \`LocalizedString\` is a branded subtype of \`string\`, it remains fully
 * backward compatible—you can pass it anywhere a \`string\` is expected.
 *
 * @example
 *   // Enforce localized strings in your components
 *   function PageTitle(props: { title: LocalizedString }) {
 *     return <h1>{props.title}</h1>
 *   }
 *
 *   // ✅ Correct: using a message function
 *   <PageTitle title={m.welcome_title()} />
 *
 *   // ❌ Type error: raw strings are not LocalizedString
 *   <PageTitle title="Welcome" />
 *
 * @example
 *   // LocalizedString is assignable to string (backward compatible)
 *   const localized: LocalizedString = m.greeting()
 *   const str: string = localized  // ✅ works fine
 *
 *   // But string is not assignable to LocalizedString
 *   const raw: LocalizedString = "Hello"  // ❌ Type error
 *
 * @example
 *   // Catches accidental string concatenation
 *   function showMessage(msg: LocalizedString) { ... }
 *
 *   showMessage(m.hello())                    // ✅
 *   showMessage("Hello " + userName)          // ❌ Type error
 *   showMessage(m.hello_user({ name: userName }))  // ✅ use params instead
 *
 * @typedef {string & { readonly __brand: 'LocalizedString' }} LocalizedString
 */
/**
 * A single markup option passed to a tag instance.
 *
 * @typedef {{
 *   name: string;
 *   value: unknown;
 * }} MessageMarkupOption
 */
/**
 * A single static markup attribute attached to a tag instance.
 *
 * @typedef {{
 *   name: string;
 *   value: string | true;
 * }} MessageMarkupAttribute
 */
/**
 * Record of markup options for a tag instance.
 *
 * @typedef {Record<string, unknown>} MessageMarkupOptions
 */
/**
 * Record of markup attributes for a tag instance.
 *
 * @typedef {Record<string, string | true>} MessageMarkupAttributes
 */
/**
 * Type-level schema for a single markup tag.
 *
 * @typedef {{
 *   options: MessageMarkupOptions;
 *   attributes: MessageMarkupAttributes;
 *   children: boolean;
 * }} MessageMarkupTag
 */
/**
 * Type-level schema for all markup tags in a message.
 *
 * @typedef {Record<string, MessageMarkupTag>} MessageMarkupSchema
 */
/**
 * Type-only metadata attached to compiled message functions.
 *
 * @template Inputs
 * @template Options
 * @template {MessageMarkupSchema} [Markup = MessageMarkupSchema]
 * @typedef {{
 *   readonly __paraglide?: {
 *     inputs: Inputs;
 *     options: Options;
 *     markup: Markup;
 *   };
 * }} MessageMetadata
 */
/**
 * A compiled, framework-neutral message part.
 *
 * @typedef {{
 *   type: "text";
 *   value: string;
 * } | {
 *   type: "markup-start";
 *   name: string;
 *   options: MessageMarkupOptions;
 *   attributes: MessageMarkupAttributes;
 * } | {
 *   type: "markup-end";
 *   name: string;
 *   options: MessageMarkupOptions;
 *   attributes: MessageMarkupAttributes;
 * } | {
 *   type: "markup-standalone";
 *   name: string;
 *   options: MessageMarkupOptions;
 *   attributes: MessageMarkupAttributes;
 * }} MessagePart
 */
/**
 * A message function is a message for a specific locale.
 *
 * @example
 *   m.hello({ name: 'world' })
 *
 * @typedef {(inputs?: Record<string, never>) => LocalizedString} MessageFunction
 */
/**
 * A message bundle function that selects the message to be returned.
 *
 * Uses `getLocale()` under the hood to determine the locale with an option.
 *
 * @template {string} T
 *
 * @example
 *   import { m } from './messages.js'
 *   m.hello({ name: 'world' }, { locale: "en" })
 *
 * @typedef {(params: Record<string, never>, options: { locale: T }) => LocalizedString} MessageBundleFunction
 */
