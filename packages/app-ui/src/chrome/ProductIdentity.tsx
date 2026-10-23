import { m, useLocale } from "../i18n/locale";

export function ProductIdentity({
  logoSrc,
  version,
  href,
}: {
  logoSrc: string;
  version: string | undefined;
  href?: string | undefined;
}) {
  const locale = useLocale();

  const content = (
    <>
      <img src={logoSrc} alt="" className="size-6" />
      <span className="font-heading text-base font-bold">TreeTime</span>
      {version !== undefined && (
        <span className="text-muted-foreground font-mono text-xs">{m.chrome_version({ version }, { locale })}</span>
      )}
    </>
  );

  return href === undefined ? (
    <div className="flex items-center gap-2 px-1">{content}</div>
  ) : (
    <a href={href} className="flex items-center gap-2 rounded-md px-1">
      {content}
    </a>
  );
}
