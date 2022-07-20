use crate::*;
use derive_more::*;

pub const VALID_CONFIGS: &[&str] = &[
    "mkl-dynamic-ilp64-iomp",
    "mkl-dynamic-ilp64-seq",
    "mkl-dynamic-lp64-iomp",
    "mkl-dynamic-lp64-seq",
    "mkl-static-ilp64-iomp",
    "mkl-static-ilp64-seq",
    "mkl-static-lp64-iomp",
    "mkl-static-lp64-seq",
];

#[derive(Debug, Clone, Copy, PartialEq, Display)]
pub enum LinkType {
    #[display(fmt = "static")]
    Static,
    #[display(fmt = "dynamic")]
    Shared,
}

#[derive(Debug, Clone, Copy, PartialEq, Display)]
pub enum Interface {
    #[display(fmt = "lp64")]
    LP64,
    #[display(fmt = "ilp64")]
    ILP64,
}

#[derive(Debug, Clone, Copy, PartialEq, Display)]
pub enum Threading {
    #[display(fmt = "iomp")]
    OpenMP,
    #[display(fmt = "seq")]
    Sequential,
}

/// Configure for linking, downloading and packaging Intel MKL
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Config {
    pub link: LinkType,
    pub index_size: Interface,
    pub parallel: Threading,
}

impl Config {
    pub fn from_str(name: &str) -> Result<Self> {
        let parts: Vec<_> = name.split("-").collect();
        if parts.len() != 4 {
            bail!("Invalid name: {}", name);
        }

        if parts[0] != "mkl" {
            bail!("Name must start with 'mkl': {}", name);
        }

        let link = match parts[1] {
            "static" => LinkType::Static,
            "dynamic" => LinkType::Shared,
            another => bail!("Invalid link spec: {}", another),
        };

        let index_size = match parts[2] {
            "lp64" => Interface::LP64,
            "ilp64" => Interface::ILP64,
            another => bail!("Invalid index spec: {}", another),
        };

        let parallel = match parts[3] {
            "iomp" => Threading::OpenMP,
            "seq" => Threading::Sequential,
            another => bail!("Invalid parallel spec: {}", another),
        };

        Ok(Config {
            link,
            index_size,
            parallel,
        })
    }

    pub fn possibles() -> Vec<Self> {
        VALID_CONFIGS
            .iter()
            .map(|name| Self::from_str(name).unwrap())
            .collect()
    }

    /// identifier used in pkg-config
    pub fn name(&self) -> String {
        format!("mkl-{}-{}-{}", self.link, self.index_size, self.parallel)
    }

    /// Common components
    ///
    /// The order must be following (or equivalent libs)
    ///
    /// mkl_intel_lp64 > mkl_intel_thread > mkl_core > iomp5
    ///
    pub fn libs(&self) -> Vec<String> {
        let mut libs = Vec::new();
        match self.index_size {
            Interface::LP64 => {
                libs.push("mkl_intel_lp64".into());
            }
            Interface::ILP64 => {
                libs.push("mkl_intel_ilp64".into());
            }
        };
        match self.parallel {
            Threading::OpenMP => {
                libs.push("mkl_intel_thread".into());
            }
            Threading::Sequential => {
                libs.push("mkl_sequential".into());
            }
        };
        libs.push("mkl_core".into());
        if matches!(self.parallel, Threading::OpenMP) {
            libs.push("iomp5".into());
        }
        libs
    }

    /// Dynamically loaded libraries, e.g. `libmkl_vml_avx2.so`
    ///
    /// - MKL seeks additional shared library **on runtime**.
    ///   This function lists these files for packaging.
    pub fn additional_libs(&self) -> Vec<String> {
        match self.link {
            LinkType::Static => Vec::new(),
            LinkType::Shared => {
                let mut libs = Vec::new();
                for prefix in &["mkl", "mkl_vml"] {
                    for suffix in &["def", "avx", "avx2", "avx512", "avx512_mic", "mc", "mc3"] {
                        libs.push(format!("{}_{}", prefix, suffix));
                    }
                }
                libs.push("mkl_rt".into());
                libs.push("mkl_vml_mc2".into());
                libs.push("mkl_vml_cmpt".into());
                libs
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_to_config() -> Result<()> {
        let cfg = Config::from_str("mkl-static-lp64-iomp")?;
        assert_eq!(
            cfg,
            Config {
                link: LinkType::Static,
                index_size: Interface::LP64,
                parallel: Threading::OpenMP
            }
        );
        Ok(())
    }

    #[test]
    fn name_to_config_to_name() -> Result<()> {
        for name in VALID_CONFIGS {
            let cfg = Config::from_str(name)?;
            assert_eq!(&cfg.name(), name);
        }
        Ok(())
    }

    #[test]
    fn invalid_names() -> Result<()> {
        assert!(Config::from_str("").is_err());
        assert!(Config::from_str("static-lp64-iomp").is_err());
        assert!(Config::from_str("mkll-static-lp64-iomp").is_err());
        assert!(Config::from_str("mkl-sttic-lp64-iomp").is_err());
        assert!(Config::from_str("mkl-static-l64-iomp").is_err());
        assert!(Config::from_str("mkl-static-lp64-omp").is_err());
        Ok(())
    }
}
