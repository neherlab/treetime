import type {
  Parsed,
  zAncestorShift,
  zAncestorState,
  zAncestralResults,
  zBranchMutations,
  zCitation,
  zClockResults,
  zCoalescentPrior,
  zDateInterval,
  zMugrationResults,
  zRecurrentSite,
  zResultNode,
  zResultTree,
  zRootToTip,
  zSettingDifference,
  zSettingsComparison,
  zSkylineSegment,
  zStateChange,
  zTimetreeEstimates,
  zTimetreeResults,
  zTreeSummary,
  zYearDate,
} from "@neherlab/app-contracts";

export type AncestorShift = Parsed<typeof zAncestorShift>;

export type AncestorState = Parsed<typeof zAncestorState>;

export type AncestralData = Parsed<typeof zAncestralResults>;

export type BranchMutations = Parsed<typeof zBranchMutations>;

export type Citation = Parsed<typeof zCitation>;

export type ClockData = Parsed<typeof zClockResults>;

export type CoalescentPrior = Parsed<typeof zCoalescentPrior>;

export type DateInterval = Parsed<typeof zDateInterval>;

export type MugrationData = Parsed<typeof zMugrationResults>;

export type RecurrentSite = Parsed<typeof zRecurrentSite>;

export type ResultNode = Parsed<typeof zResultNode>;

export type ResultTree = Parsed<typeof zResultTree>;

export type RootToTip = Parsed<typeof zRootToTip>;

export type SettingDifference = Parsed<typeof zSettingDifference>;

export type SettingsComparison = Parsed<typeof zSettingsComparison>;

export type SkylineSegment = Parsed<typeof zSkylineSegment>;

export type StateChange = Parsed<typeof zStateChange>;

export type TimetreeData = Parsed<typeof zTimetreeResults>;

export type TimetreeEstimates = Parsed<typeof zTimetreeEstimates>;

export type TreeSummary = Parsed<typeof zTreeSummary>;

export type YearDate = Parsed<typeof zYearDate>;
