# PharmFlux R source package

The R package uses the shared native Rust engine. Compile a model once with
`compile_model()`, then simulate, fit, scan, or run Morris screening using
versioned requests. Results retain units, observation sides, and execution
identity. Check fitting status before interpreting parameter estimates.

Start with the [R user guide](../../docs/R-USER-GUIDE.md), or download its
[R Markdown tutorial](../../docs/R-USER-GUIDE.Rmd). Run examples against the
same source checkout used to build the package.

## Install from the public workspace

R, a compatible C compiler, Rust 1.98 or later, and Cargo are required.
From the repository root on macOS or Linux:

```sh
mkdir -p "$HOME/R/pharmflux-library"
export R_LIBS_USER="$HOME/R/pharmflux-library"
export R_HOME="$(R RHOME)"
cargo fetch --locked
R CMD INSTALL --library="$R_LIBS_USER" bindings/r
Rscript -e 'library(pharmflux); print(packageVersion("pharmflux"))'
```

The package builds its native library using the surrounding public Rust
workspace. Keep the entire checkout together. Registry dependencies must be
cached for its locked offline Cargo build. Restart R after replacing a loaded
native library.

A standalone source archive can instead be generated with its own workspace:

```sh
python3 scripts/prepare-r-source.py /tmp/pharmflux-r-source
R CMD INSTALL --library="$R_LIBS_USER" /tmp/pharmflux-r-source/pharmflux_0.1.0.tar.gz
```

The source archive still needs the cached Cargo registry dependencies. These
instructions do not assume a CRAN or R-universe release. The documentation
examples are checked on macOS ARM64; CI checks the Linux source installation.
Windows installation requires separate verification.

## Data and estimation boundaries

`nonmem_requests()` handles an explicit subset of NONMEM-style records with
units, dose compartments, and observation sides supplied by the caller.
Unsupported records are rejected. `fit_data_requests()` prepares individual
or pooled Gaussian fits. `population_fit_data_requests()` prepares the native
FOCEI or SAEM request and retains source-row mappings. These helpers are not
a general NONMEM importer: they do not infer covariates, units, endpoints, or
missing-data policy.

The package includes installed R help for each exported function. See the
user guide for complete examples, fitting limits, and error handling. The
engine and R wrappers use Apache-2.0; dependency licenses remain separate.
