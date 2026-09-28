use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

const CITATION_TEXT: &str = "Sagulenko P, Puller V, Neher RA. TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution 4 (2018), vex042.";

const CITATION_DOI: &str = "10.1093/ve/vex042";

const DOI_RESOLVER: &str = "https://doi.org/";

/// The publication to cite for TreeTime.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct Citation {
  /// Reference in text form.
  pub text: String,
  /// DOI of the publication, for example `10.1093/ve/vex042`.
  pub doi: String,
  /// Link to the publication.
  pub url: String,
}

pub fn citation() -> Citation {
  Citation {
    text: CITATION_TEXT.to_owned(),
    doi: CITATION_DOI.to_owned(),
    url: format!("{DOI_RESOLVER}{CITATION_DOI}"),
  }
}
