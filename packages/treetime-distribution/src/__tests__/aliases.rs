use crate::Distribution;
use crate::policy::{NegLog, Plain};

pub(crate) type DistributionPlain = Distribution<Plain>;
pub(crate) type DistributionNegLog = Distribution<NegLog>;
