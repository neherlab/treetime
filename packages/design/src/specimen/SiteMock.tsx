import { Avatar, AvatarFallback } from "../ui/avatar";
import { Wordmark } from "../ui/wordmark";

const NAV = ["Research", "Publications", "Software", "Talks", "Teaching", "Team"] as const;

const AREAS = [
  { title: "Evolution of bacteria", topics: ["horizontal transfer", "pangenomes", "drug resistance"] },
  { title: "HIV evolution", topics: ["evolution within patients", "mutation and recombination", "latent reservoir"] },
  { title: "Molecular epidemiology", topics: ["real-time tracking", "prediction", "antigenic evolution"] },
  { title: "Theory", topics: ["population genetics", "rapid adaptation", "time tree inference"] },
] as const;

const TOOLS = [
  { name: "TreeTime", use: "Time-scaled phylogenies" },
  { name: "Nextclade", use: "Clade assignment and quality control" },
  { name: "Nextstrain", use: "Real-time tracking of pathogen evolution" },
  { name: "Pangraph", use: "Pangenome graphs of bacteria" },
] as const;

const PUBLICATIONS = [
  {
    title: "Sample title: epistasis shapes the fitness landscape of a virus",
    authors: "Author A and Author B",
    venue: "Sample Journal 12, 345",
    date: "17 September 2026",
  },
  {
    title: "Sample title: antibody landscapes inform the choice of vaccine strains",
    authors: "Author C, Author D, and others",
    venue: "Sample Journal 9, 46",
    date: "26 July 2026",
  },
] as const;

const TEAM = [
  { initials: "PI", role: "Principal investigator" },
  { initials: "SS", role: "Staff scientist" },
  { initials: "PD", role: "Postdoc" },
  { initials: "SD", role: "Software developer" },
] as const;

export function SiteMock() {
  return (
    <div className="bg-card overflow-hidden rounded-xl border">
      <header className="flex h-14 items-center gap-6 border-b px-6">
        <Wordmark name="neherlab" place="biozentrum" />
        <nav aria-label="Lab site" className="ml-auto hidden gap-5 text-sm @3xl:flex">
          {NAV.map((item) => (
            <a key={item} href="#site" className="hover:text-foreground text-muted-foreground first:text-foreground">
              {item}
            </a>
          ))}
        </nav>
      </header>
      <div className="mx-auto flex max-w-5xl flex-col gap-12 px-6 py-12">
        <div className="flex flex-col gap-5 border-b pb-10">
          <h3 className="text-brand font-serif text-4xl font-normal tracking-tight text-balance @3xl:text-6xl">
            Pathogen evolution, genomics, and biophysics
          </h3>
          <p className="max-w-[60ch] text-lg leading-relaxed">
            We are a research group at the{" "}
            <a href="#site" className="text-primary underline underline-offset-3">
              Biozentrum of the University of Basel
            </a>
            . We study how pathogens such as HIV, influenza viruses, and bacteria evolve, and we build open tools that
            follow their spread.
          </p>
        </div>
        <ResearchAreas />
        <Publications />
        <Team />
      </div>
      <footer className="text-muted-foreground flex justify-between border-t px-6 py-4 text-xs">
        <span>Biozentrum, University of Basel</span>
        <a href="#site" className="underline">
          Imprint
        </a>
      </footer>
    </div>
  );
}

function ResearchAreas() {
  return (
    <div className="grid gap-10 @4xl:grid-cols-[1fr_16rem]">
      <div className="grid gap-8 @2xl:grid-cols-2">
        {AREAS.map((area) => (
          <div key={area.title} className="flex flex-col gap-2">
            <h4 className="font-serif text-2xl font-normal">
              <a href="#site" className="text-primary hover:underline">
                {area.title}
              </a>
            </h4>
            <ul className="list-disc pl-5 text-sm leading-relaxed">
              {area.topics.map((topic) => (
                <li key={topic}>{topic}</li>
              ))}
            </ul>
          </div>
        ))}
      </div>
      <div className="flex flex-col gap-3">
        <h4 className="font-serif text-2xl font-normal">Tools</h4>
        <ul className="flex flex-col gap-2.5">
          {TOOLS.map((tool) => (
            <li key={tool.name} className="flex flex-col">
              <a href="#site" className="text-primary hover:underline">
                {tool.name}
              </a>
              <span className="text-muted-foreground text-xs">{tool.use}</span>
            </li>
          ))}
        </ul>
      </div>
    </div>
  );
}

function Publications() {
  return (
    <div className="flex flex-col gap-2">
      <div className="flex items-baseline justify-between gap-4">
        <h4 className="font-serif text-4xl font-normal tracking-tight">Recent publications</h4>
        <a href="#site" className="text-muted-foreground hover:text-primary font-serif text-xl">
          All publications
        </a>
      </div>
      <ul>
        {PUBLICATIONS.map((publication) => (
          <li
            key={publication.title}
            className="border-input grid grid-cols-[1fr_auto] gap-x-6 gap-y-0.5 border-b border-dotted py-3"
          >
            <a href="#site" className="text-primary hover:underline">
              {publication.title}
            </a>
            <span className="text-muted-foreground col-start-1 text-sm">{publication.authors}</span>
            <span className="text-muted-foreground col-start-1 text-sm">{publication.venue}</span>
            <span className="text-muted-foreground col-start-2 row-span-3 row-start-1 self-end text-sm">
              {publication.date}
            </span>
          </li>
        ))}
      </ul>
    </div>
  );
}

function Team() {
  return (
    <div className="flex flex-col gap-5">
      <h4 className="font-serif text-4xl font-normal tracking-tight">Team</h4>
      <ul className="grid grid-cols-2 gap-6 @3xl:grid-cols-4">
        {TEAM.map((person) => (
          <li key={person.initials} className="flex flex-col items-start gap-2">
            <Avatar size="lg">
              <AvatarFallback>{person.initials}</AvatarFallback>
            </Avatar>
            <span className="text-primary text-sm font-medium">Name Surname</span>
            <span className="text-muted-foreground text-xs">{person.role}</span>
          </li>
        ))}
      </ul>
    </div>
  );
}
