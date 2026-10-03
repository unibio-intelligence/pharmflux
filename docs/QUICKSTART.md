# PharmFlux user guides

Use the executable guide for your analysis language. Each starts with a
self-contained one-compartment simulation and then covers oral models, dosing,
result tables, numerical settings, scans, fitting, population data,
reproducibility, and troubleshooting. All examples use synthetic data and the
source-based 0.1.0 bindings.

| Language | Executable guide | Rendered guide | Plain guide | Script |
| --- | --- | --- | --- | --- |
| Python | [Jupyter notebook](PYTHON-USER-GUIDE.ipynb) | [HTML with outputs](PYTHON-USER-GUIDE.html) | [Markdown](PYTHON-USER-GUIDE.md) | [Python script](examples/python_quickstart.py) |
| R | [R Markdown](R-USER-GUIDE.Rmd) | [HTML with outputs](R-USER-GUIDE.html) | [Markdown](R-USER-GUIDE.md) | [R script](examples/r_quickstart.R) |

## Run the executable guides

Install PharmFlux using the corresponding guide before execution. For Python,
also install notebook dependencies in the same environment:

```sh
python -m pip install jupyterlab ipykernel pandas matplotlib nbformat nbconvert
jupyter lab docs/PYTHON-USER-GUIDE.ipynb
```

Select the kernel that has PharmFlux installed, restart it, and run all cells
in order. The notebook retains checked outputs and an inline plot. New analysis
files go into `pharmflux-python-tutorial-output` beneath the kernel's starting
directory. Installation commands are reference text rather than executed cells.

For R, install rendering dependencies, then render from the checkout root:

```r
install.packages(c("rmarkdown", "knitr"))
# Ensure the installed PharmFlux library is on .libPaths() before rendering.
rmarkdown::render("docs/R-USER-GUIDE.Rmd")
```

Pandoc is also required; RStudio commonly supplies it, or use a separate Pandoc
installation. The Rmd executes its analysis chunks in order and displays plots,
tables and fitting results in HTML. Installation chunks are marked `eval=FALSE`.
Analysis files go into `docs/pharmflux-r-tutorial-output` when rendering the
file at its supplied location.

## Execution checks

Both formats were executed using installed local PharmFlux bindings. The
Jupyter run completed all 18 code cells with zero error outputs, including
the pandas/matplotlib plot. The Rmd rendered successfully with 17 analysis
example chunks plus its setup chunk, and an embedded plot. The checks use fresh installed bindings built from the source being published.

The checks cover analytic simulation, event sides, oral dose handling, scans,
Morris screening, individual parameter recovery, a focused population example,
and output saving. The short SAEM demonstration intentionally uses a small
iteration budget and reports its status; it is not a recommended estimation
configuration. These checks do not rerun the installation instructions or
establish clean source installation on every platform.

## Maintain the guides

The Markdown guides supply the narrative and code for both executable formats.
After editing them, regenerate the notebook and Rmd from the checkout root:

```sh
python scripts/build_user_guide_formats.py
```

Generation creates an unexecuted notebook and an Rmd. Rerun the notebook and
render the Rmd before replacing checked HTML or updating the execution record.
This avoids maintaining separate copies of the tutorial content by hand.

The standalone scripts remain useful for a dependency-light workflow:

```sh
python docs/examples/python_quickstart.py
Rscript docs/examples/r_quickstart.R
```

They write a trajectory and `pharmflux-analysis` directory to the current
working directory. The R script also saves a plot. The Python script requires
only PharmFlux and the standard library.

For related API details, see [Python binding notes](PYTHON.md),
[R binding notes](../bindings/r/README.md), and [SAEM scope](SAEM.md).
